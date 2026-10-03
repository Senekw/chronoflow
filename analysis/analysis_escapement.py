"""Escapement + balance dynamics for ChronoFlow.

Balance:  I th'' = -k th - c th' + tau_esc(th, th', direction)
Two escapement models, both piecewise-constant torque in balance-angle windows:

(A) Conventional pin-pallet LEVER escapement (baseline).  Per beat the escape wheel
    advances half a tooth pitch and releases E_beat = T_e * pi / Z.  The balance receives
    E_imp = IMPULSE_FRACTION*E_beat over the lift arc (after unlocking) and pays
    E_unl = UNLOCK_FRACTION*E_beat over the unlocking arc.  Both scale with T_e.

(B) Constant-force SPRING-ARM escapement (ChronoFlow).  A spring-loaded impulse arm,
    cocked by the escape wheel and held by a latch, is released by the balance at
    -ARM_RELEASE and pushes the impulse pin over ARM_IMPULSE delivering a fixed E_IMPULSE;
    the balance pays a fixed E_UNLATCH to trip the latch.  Train torque only re-cocks the
    arm; it must exceed T_min = (E_IMPULSE/LIFT_EFF)/(pi/Z) or the escapement stalls.

Period is measured over steady-state zero crossings; the amplitude is the mean peak angle.
"""
import numpy as np
import design_params as dp


def simulate(I, k, Q, model, T_e=None, t_end=120.0, dt=2e-4, theta0=0.5, trace=False):
    w0 = np.sqrt(k / I)
    c = I * w0 / Q
    Z = dp.Z_ESCAPE
    if model == "lever":
        E_beat = T_e * np.pi / Z
        E_imp = dp.IMPULSE_FRACTION * E_beat
        E_unl = dp.UNLOCK_FRACTION * E_beat
        th_u = np.deg2rad(dp.UNLOCK_ANGLE_DEG)
        th_i = np.deg2rad(dp.LIFT_ANGLE_DEG)
        half = 0.5 * (th_u + th_i)
        entry = -half                     # progress coordinate u = s*theta
        unl_end = entry + th_u
        imp_end = half
    elif model == "arm":
        T_min = (dp.E_IMPULSE / dp.LIFT_EFF) / (np.pi / Z)
        if T_e is not None and T_e < T_min:
            return dict(status="stall", reason=f"T_e {T_e:.2e} < T_min {T_min:.2e}", T_min=T_min)
        E_imp = dp.E_IMPULSE
        E_unl = dp.E_UNLATCH
        th_r = np.deg2rad(dp.ARM_RELEASE_DEG)
        th_i = np.deg2rad(dp.ARM_IMPULSE_DEG)
        th_u = np.deg2rad(2.0)            # latch tripping arc
        entry = -th_r - th_u
        unl_end = -th_r
        imp_end = -th_r + th_i
    else:
        raise ValueError(model)
    tau_i = E_imp / (imp_end - unl_end)
    tau_u = E_unl / (unl_end - entry)

    n = int(t_end / dt)
    th, om = theta0, 0.0
    engaged = False          # pin inside fork/arm window
    stalled = False
    peaks, zc_times = [], []
    tr_t, tr_th = [], []
    last_sign = np.sign(th)
    prev_om = om
    t = 0.0

    def tau_esc(th, om, engaged):
        if not engaged:
            return 0.0
        s = 1.0 if om > 0 else -1.0
        u = s * th
        if u < entry or u > imp_end:
            return 0.0
        if u < unl_end:
            return -s * tau_u
        return s * tau_i

    def f(th, om, engaged):
        return om, (-k * th - c * om + tau_esc(th, om, engaged)) / I

    for i in range(n):
        # engagement bookkeeping: pin enters window when |th| crosses 'entry' moving inward
        s = 1.0 if om > 0 else -1.0
        u = s * th
        if not engaged and entry <= u <= imp_end and (s * prev_om > 0):
            engaged = True
        if engaged and u > imp_end:
            engaged = False
        # RK4
        k1 = f(th, om, engaged)
        k2 = f(th + 0.5*dt*k1[0], om + 0.5*dt*k1[1], engaged)
        k3 = f(th + 0.5*dt*k2[0], om + 0.5*dt*k2[1], engaged)
        k4 = f(th + dt*k3[0], om + dt*k3[1], engaged)
        th_new = th + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        om_new = om + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        # reversal while the pin is still inside the fork/arm window -> the escapement
        # cannot complete unlocking + impulse (wheel never re-locks) -> fault/stall
        if engaged and (om_new * om < 0):
            stalled = True
            break
        # peaks
        if om_new * om < 0:
            peaks.append((t, abs(th_new)))
        # positive-going zero crossings
        if th < 0 <= th_new:
            zc_times.append(t + dt * (-th) / (th_new - th))
        prev_om = om
        th, om = th_new, om_new
        t += dt
        if trace and i % 25 == 0:
            tr_t.append(t); tr_th.append(th)
        if abs(th) > dp.THETA_MAX * 1.5:
            out = dict(status="overbank", amplitude=abs(th))
            if trace: out["t"] = np.array(tr_t); out["theta"] = np.array(tr_th)
            return out
    if stalled:
        return dict(status="stall", reason="reversal inside engagement window", t=t)
    zc = np.array(zc_times)
    pk = np.array(peaks)
    sel = zc > 0.6 * t_end
    if sel.sum() < 5:
        return dict(status="nosteady")
    periods = np.diff(zc[sel])
    amp = pk[pk[:, 0] > 0.6 * t_end, 1].mean()
    T = periods.mean()
    out = dict(status="ok", period=T, freq=1 / T, period_jitter=periods.std() / T,
               amplitude=amp, free_period=2*np.pi/w0,
               rate_error=(T - 2*np.pi/w0) / (2*np.pi/w0))
    if trace:
        out["t"] = np.array(tr_t); out["theta"] = np.array(tr_th)
    return out


def escape_wheel_torque(F_excess, flow_mL_h, area_mm2=dp.SYR_AREA):
    """Torque available at the escape wheel from the carriage spring excess force,
    through rack + fixed train + cassette (efficiencies), at oscillator design frequency."""
    v = flow_mL_h * 1000 / area_mm2 / 3600 * 1e-3      # m/s
    P = F_excess * v
    eta = dp.RACK_EFF * dp.GEAR_EFF ** (len(dp.STAGES) + 2)
    w_e = 2 * np.pi * (2 * dp.F_TARGET / dp.BEATS_PER_REV)
    return P * eta / w_e, P, eta


if __name__ == "__main__":
    import analysis_oscillator as ao
    k = ao.k_total(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0, dp.BAND_T0)
    m, I = ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG,
                              dp.N_TRIM, dp.M_TRIM, dp.R_TRIM, dp.HUB_M, dp.HUB_R)
    for Te in (2e-6, 5e-6, 2e-5, 1e-4):
        print("lever", Te, simulate(I, k, 150, "lever", T_e=Te, t_end=80))
    print("arm", simulate(I, k, 150, "arm", T_e=2e-5, t_end=80))
