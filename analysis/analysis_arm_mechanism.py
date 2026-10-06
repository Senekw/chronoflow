"""Mechanism-level analysis of the ChronoFlow constant-force spring-arm escapement (rev 4, review round 5).

The beat-level simulation (analysis_escapement.simulate, model "arm") assumes a fixed impulse energy E_IMPULSE and a
fixed unlatch energy E_UNLATCH.  This module derives both from the mechanism and quantifies what the 24x train-torque
range does to the escape wheel itself:

  1. Impulse energy and its tolerance  -- preloaded arm spring unwinding over the impulse arc; latch-position tolerance
  2. Unlatch work                      -- latch normal load x friction x engagement depth, mu = 0.1 ... 0.4
  3. Wheel force vs arm force          -- shows that the wheel must be held by a separate detent during impulse
  4. Re-cocking dynamics (ODE)         -- time for the wheel to run half a pitch and re-cock the arm vs the one-period budget
  5. Lock impact                       -- surplus energy per beat, impact velocity, peak force and Hertz stress at the detent,
                                          compliant (printed flexure) vs rigid stop; impact count over a 20 h infusion
  6. Sensitivity of the rate           -- simulate() with E_IMPULSE +/-10 % and E_UNLATCH over the friction range

Outputs: results/arm_mechanism.json, figures/fig7_arm_mechanism.png
"""
import os, json, math, importlib
import numpy as np
import matplotlib.pyplot as plt
import design_params as dp, analysis_escapement as ae
importlib.reload(dp); importlib.reload(ae)
try:
    apply_figure_style
except NameError:
    import os as _os, sys as _sys; _sys.path.append(_os.path.dirname(_os.path.abspath(__file__)) if "__file__" in globals() else _os.getcwd())
    from figstyle import apply_figure_style, panel_letter, META_GREY

os.makedirs("results", exist_ok=True); os.makedirs("figures", exist_ok=True)
res = {}

# ------------------------------------------------------------------ geometry (CAD representative embodiment, cad_model.spring_arm_escapement)
Z = dp.Z_ESCAPE; HALF_PITCH = math.pi / Z                 # rad of wheel per beat
R_W = 12.0e-3                                            # m, escape-wheel tip radius
R_P = 5.0e-3                                             # m, arm pivot -> pallet pin
R_T = 59.0e-3                                            # m, arm pivot -> impulse tip
R_PIN = dp.HUB_R + 1.2e-3                                # m, impulse-pin radius on the balance hub (6.2 mm)
R_LATCH = 16.0e-3                                        # m, arm pivot -> latch
D_LATCH = 0.3e-3                                         # m, latch engagement depth (sliding distance to release)
TH_IMP = math.radians(dp.ARM_IMPULSE_DEG)                # balance arc of impulse
S_IMP = R_PIN * TH_IMP                                   # tip travel during impulse (1.62 mm)
DTH_ARM = S_IMP / R_T                                    # arm rotation during impulse (27.5 mrad)
PHI_COCK = 0.10                                          # rad of wheel rotation over which the tooth back-slope re-cocks the arm
PHI_FREE = HALF_PITCH - PHI_COCK                         # free-running arc before the tooth reaches the pallet
F_BEAT = 2 * dp.F_TARGET                                 # beats per second
PERIOD = 1.0 / dp.F_TARGET
T_IMPULSE = TH_IMP / (2 * math.pi * dp.F_TARGET * 0.6)   # ~ time the pin spends in the impulse arc at amplitude 0.6 rad (order of magnitude)
# rev 5: the two arms alternate, so each ARM has a full period to be re-cocked, but the single escape WHEEL must run its
# half pitch (re-cock + lock) between successive beats, i.e. within half a period minus the impulse time.  Event timeline per
# beat: unlatch -> impulse (T_IMPULSE) -> arm end-of-stroke releases the detent -> wheel free arc -> cocking arc -> lock.
BUDGET = PERIOD / 2 - T_IMPULSE
BUDGET_FULL_PERIOD = PERIOD - T_IMPULSE                  # the (incorrect) rev-4 budget, kept for the record

# ------------------------------------------------------------------ 1. impulse energy: FORWARD model from a specified spring
# rev 5 (review round 6): the energy is no longer taken from design_params.E_IMPULSE; it is the work a specified torsion
# spring (k, preload) does on the impulse pin as the arm follows the pin over the impulse arc, integrated through the
# contact geometry (tip on a pin that moves on a circle of radius R_PIN about the balance axis), minus tip friction.
K_F = 1.0e-4               # N m / rad, music-wire torsion spring on the arm pivot (specified)
TH_PRE = math.radians(26.3)  # rad, preload angle at the start of the impulse (set by the latch position); chosen so the forward model gives 1.2 uJ
MU_TIP = 0.15              # steel tip on steel pin, lightly oiled
def impulse_energy(k=K_F, th_pre=TH_PRE, mu=MU_TIP, n=400):
    """Work delivered to the balance over the impulse arc (balance angle -ARM_RELEASE .. -ARM_RELEASE+ARM_IMPULSE).
    The arm axis passes through the balance centre, so the tip pushes the pin along the pin's tangent at the centre
    crossing; off-centre the tip force makes angle phi with the pin's tangent (cos phi work factor) and the sliding
    component mu*N*sin(phi) is lost.  Arm rotation d(theta) = R_PIN cos(phi) d(phi) / R_T."""
    phi0 = -math.radians(dp.ARM_RELEASE_DEG); phi1 = phi0 + TH_IMP
    phis = np.linspace(phi0, phi1, n + 1)
    th = 0.0; E = 0.0
    for a, b in zip(phis[:-1], phis[1:]):
        dphi = b - a; phim = 0.5 * (a + b)
        dth = R_PIN * math.cos(phim) * dphi / R_T                      # arm rotation for this pin step
        tau = k * (th_pre - th)                                        # arm torque at the current unwinding
        F_tip = tau / R_T                                              # force at the tip, normal to the arm
        dE = F_tip * math.cos(phim) * (R_PIN * dphi) - mu * F_tip * abs(math.sin(phim)) * (R_PIN * dphi)
        E += dE; th += dth
    return E, th
E_imp, DTH_ARM_FWD = impulse_energy()
res["impulse"] = dict(spring_k_Nm_per_rad=K_F, preload_deg=math.degrees(TH_PRE), mu_tip=MU_TIP,
                      E_impulse_forward_J=E_imp, E_impulse_assumed_in_beat_model_J=dp.E_IMPULSE,
                      arm_rotation_mrad=DTH_ARM_FWD*1e3, tip_travel_mm=S_IMP*1e3,
                      E_vs_preload={f"{d:+d} deg": impulse_energy(th_pre=TH_PRE + math.radians(d))[0] for d in (-3, -2, -1, 0, 1, 2, 3)},
                      E_vs_mu={f"{m:.2f}": impulse_energy(mu=m)[0] for m in (0.0, 0.1, 0.15, 0.2, 0.3)},
                      force_droop_over_impulse_pct=100 * DTH_ARM_FWD / TH_PRE)
# derived quantities used downstream
TAU_ARM = K_F * (TH_PRE - DTH_ARM_FWD / 2)                             # mean arm torque over the stroke
DTH_ARM = DTH_ARM_FWD
F_tip = TAU_ARM / R_T; F_PALLET = TAU_ARM / R_P
LATCH_TOL = 0.1e-3; LATCH_SCATTER = 5e-6
res["impulse"].update(dict(mean_arm_torque_Nm=TAU_ARM, tip_force_mN=F_tip*1e3, pallet_force_mN=F_PALLET*1e3,
                           energy_bias_from_latch_tolerance_pct=100 * (LATCH_TOL / R_LATCH) / TH_PRE,
                           beat_to_beat_energy_scatter_pct=100 * (LATCH_SCATTER / R_LATCH) / TH_PRE))
# ---- why not the printed flexure: k = 3EI/L^3 * r^2 for the 0.8 x 2.5 x 14 mm PETG cantilever acting at 8 mm
K_FLEX = 3 * 2.0e9 * (2.5e-3 * 0.8e-3**3 / 12) / (14e-3)**3 * (8e-3)**2       # = 0.0149 N m/rad (150x stiffer than K_F)
# minimum preload that keeps the force positive to the end of the stroke is th_pre = DTH_ARM (100 % droop); the energy
# it then delivers is k*DTH^2/2; to deliver the design 1.2 uJ with positive force the required preload is below DTH -> impossible
E_flex_min_pos = 0.5 * K_FLEX * DTH_ARM**2
th_pre_flex_for_target = dp.E_IMPULSE / (K_FLEX * DTH_ARM) + DTH_ARM / 2
res["impulse"]["printed_flexure"] = dict(k_Nm_per_rad=K_FLEX, stiffness_ratio_vs_design=K_FLEX / K_F,
    preload_for_1p2uJ_rad=th_pre_flex_for_target, arm_stroke_rad=DTH_ARM,
    feasible_with_positive_force=bool(th_pre_flex_for_target >= DTH_ARM),
    energy_at_100pct_droop_uJ=E_flex_min_pos*1e6,
    note="the flexure is 150x too STIFF, not too soft: at the design energy its force would pass through zero before the "
         "end of the 27 mrad stroke (preload 16.7 mrad < stroke); a soft wire spring with a large preload is needed "
         "for a near-constant force (rev 4 text had the preload angle 1000x wrong)")

# ------------------------------------------------------------------ 2. unlatch work
N_LATCH = TAU_ARM / R_LATCH
E_latch_spring = 0.02e-6                                 # J, latch return-spring work over D_LATCH (light spring)
MU = np.array([0.1, 0.2, 0.3, 0.4])
E_unl = MU * N_LATCH * D_LATCH + E_latch_spring
res["unlatch"] = dict(latch_normal_load_mN=N_LATCH*1e3, engagement_mm=D_LATCH*1e3, mu=MU.tolist(), E_unlatch_uJ=(E_unl*1e6).tolist(),
                      E_unlatch_assumed_uJ=dp.E_UNLATCH*1e6, fraction_of_impulse_pct=(100*E_unl/E_imp).tolist())

# ------------------------------------------------------------------ 3. wheel force on the pallet vs arm force -> detent required
T_e_map = np.array([19.5e-6, 476e-6])                    # operating-map extremes of escape-wheel torque
F_wheel = T_e_map / R_W
res["lock"] = dict(wheel_force_at_pallet_mN=(F_wheel*1e3).tolist(), arm_pallet_force_mN=F_PALLET*1e3,
                   torque_above_which_wheel_overpowers_arm_uNm=F_PALLET * R_W * 1e6,
                   conclusion="the wheel force exceeds the arm force above %.0f uN m, i.e. over most of the operating map: the wheel must be "
                              "locked by a separate detent during impulse and released by the arm at the end of its stroke" % (F_PALLET*R_W*1e6))

# ------------------------------------------------------------------ 4. re-cocking dynamics
def J_eff(cassette):
    """Inertia seen at the escape wheel: wheel (brass) + g4 + upstream gears reflected by their speed ratios (solid-disk estimates)."""
    z1, z2, z3, z4 = dp.CASSETTES[cassette]
    def disk(z, rho=dp.RHO_PETG, face=dp.GEAR_FACE):
        r = 1e-3 * z / 2 + 1e-3; m = rho * math.pi * r**2 * (face * 1e-3); return 0.5 * m * r**2   # face in mm
    J = 5.70e-7                                          # escape wheel, brass, from the CAD solid
    J += disk(z4)
    n = z4 / z3                                           # S6 speed / S7 speed
    J += (disk(z3) + disk(z2)) * n**2
    n *= z2 / z1
    J += (disk(z1) + disk(12)) * n**2                    # S5: g1 + P4
    n *= 12 / 30
    J += (disk(30) + disk(12)) * n**2                    # S4: W4 + P3 (and so on: negligible beyond)
    return J

def recock(T_e, J, c_visc=2e-7):
    """Wheel runs HALF_PITCH from rest: free arc, then cocking arc with load torque T_COCK; returns time and lock speed."""
    T_COCK = (TAU_ARM * DTH_ARM / dp.LIFT_EFF) / PHI_COCK   # mean torque the tooth back-slope needs to re-cock the arm (incl. losses)
    phi, om, t, dt = 0.0, 0.0, 0.0, 1e-4
    while phi < HALF_PITCH and t < 5.0:
        load = T_COCK if phi >= PHI_FREE else 0.0
        a = (T_e - load - c_visc * om) / J
        om_new = om + a * dt
        if om_new < 0:                                   # cannot push through the cocking arc
            return dict(time_s=float("inf"), omega_lock=0.0, stalled=True)
        phi += 0.5 * (om + om_new) * dt; om = om_new; t += dt
    return dict(time_s=t, omega_lock=om, stalled=t >= 5.0)

T_grid = np.logspace(math.log10(10e-6), math.log10(500e-6), 40)
# per-cassette torque range from the operating map (torque at the escape wheel scales with flow x excess force)
try:
    import pandas as pd
    _om = pd.read_csv("results/operating_map.csv").groupby("cassette")["T_e_uNm"].agg(["min", "max"])
    T_RANGE = {c: (float(_om.loc[c, "min"])*1e-6, float(_om.loc[c, "max"])*1e-6) for c in _om.index}
except Exception:
    T_RANGE = {"1 mL/h": (19.47e-6, 47.56e-6), "2 mL/h": (38.94e-6, 95.12e-6), "5 mL/h": (97.36e-6, 237.8e-6), "10 mL/h": (194.7e-6, 475.6e-6)}
dyn = {}
for cas in dp.CASSETTES:
    J = J_eff(cas)
    rows = [recock(T, J) for T in T_grid]
    lo, hi = T_RANGE[cas]
    dyn[cas] = dict(J_eff_kgm2=J, T_uNm=(T_grid*1e6).tolist(), t_recock_s=[r["time_s"] for r in rows],
                    omega_lock=[r["omega_lock"] for r in rows], T_range_uNm=[lo*1e6, hi*1e6],
                    t_recock_at_min_torque_s=recock(lo, J)["time_s"], t_recock_at_max_torque_s=recock(hi, J)["time_s"])
    dyn[cas]["time_margin_at_min_torque"] = BUDGET / dyn[cas]["t_recock_at_min_torque_s"]
res["recock"] = dict(period_s=PERIOD, beat_interval_s=PERIOD/2, impulse_time_s=T_IMPULSE, budget_s=BUDGET, budget_rev4_full_period_s=BUDGET_FULL_PERIOD,
                     timeline=["balance pin trips latch (unlatch work)", "arm impulse over 15 deg (T_IMPULSE)", "arm tail releases the wheel detent via the release star",
                               "wheel free arc PHI_FREE", "wheel cocking arc PHI_COCK re-latches the other arm", "detent re-locks the wheel (impact)", "next beat: other arm"],
                     phi_free_rad=PHI_FREE, phi_cock_rad=PHI_COCK,
                     T_cock_uNm=(TAU_ARM * DTH_ARM / dp.LIFT_EFF) / PHI_COCK * 1e6, cassettes=dyn)
res["recock"]["worst_case"] = min(((cas, d["t_recock_at_min_torque_s"], d["time_margin_at_min_torque"]) for cas, d in dyn.items()), key=lambda x: x[2])

# ------------------------------------------------------------------ 4b. detent release kinematics (rev 5, from the CAD geometry)
TAIL_A, TAIL_B = 16.0e-3, 14.0e-3                      # arm release tails (A ends 1 mm beyond the wheel axis, B 1 mm short)
PEG_R = math.hypot(1.0e-3, 9.0e-3 + 1.95e-3)           # star peg radius from the escape axis
DET_PIN_R, DET_ARM, DET_TAIL = 11.9e-3, 12.0e-3, 4.0e-3
T_vec = np.array([DET_PIN_R - DET_TAIL, DET_ARM])      # L-tail end in (radial, tangential) coordinates at the pin
T_r = float(np.hypot(*T_vec))
cos_fac = T_vec[0] / T_r                               # fraction of the star-peg motion that moves the tail along its own motion direction
rel = {}
for name, Lt in (("arm A", TAIL_A), ("arm B", TAIL_B)):
    tail_move = Lt * DTH_ARM
    star_rot = tail_move / PEG_R
    tail_push = star_rot * T_r * cos_fac
    pin_withdraw = tail_push * DET_ARM / DET_TAIL
    rel[name] = dict(tail_motion_mm=tail_move*1e3, star_rotation_mrad=star_rot*1e3, detent_tail_push_mm=tail_push*1e3, pin_withdrawal_mm=pin_withdraw*1e3)
engage = (R_W + 0.5e-3) - DET_PIN_R                    # pin must clear the tip circle by its own radius
res["detent_release"] = dict(per_arm=rel, required_withdrawal_mm=engage*1e3,
                             margin={k: v["pin_withdrawal_mm"] / (engage*1e3) for k, v in rel.items()},
                             note="first-order kinematics of the drawn linkage; the small margin is a bench-tuning item (tail length, pin radius)")

# ------------------------------------------------------------------ 5. lock impact at the detent
K_STOP = {"rigid steel pin on PETG (1e6 N/m)": 1.0e6, "printed flexure detent 0.8x4x8 mm (2000 N/m)": 3 * 2.0e9 * (4e-3 * 0.8e-3**3 / 12) / (8e-3)**3}
E_star = 1.0 / ((1 - 0.3**2) / 200e9 + (1 - 0.4**2) / 2.0e9)     # steel on PETG
L_C, R_C = 2.5e-3, 0.5e-3                                       # pin contact length and radius
def impact(T_e, J):
    E_surplus = T_e * HALF_PITCH - TAU_ARM * DTH_ARM / dp.LIFT_EFF   # kinetic energy arriving at the lock (viscous loss neglected)
    E_surplus = max(E_surplus, 0.0)
    om = math.sqrt(2 * E_surplus / J)
    out = dict(E_surplus_uJ=E_surplus*1e6, omega_lock=om, rim_speed_mm_s=om*R_W*1e3)
    for name, k in K_STOP.items():
        F = math.sqrt(2 * E_surplus * k); x = math.sqrt(2 * E_surplus / k) if k else 0
        p = math.sqrt(F * E_star / (math.pi * L_C * R_C)) if F > 0 else 0.0
        out[name] = dict(peak_force_N=F, deflection_mm=x*1e3, hertz_stress_MPa=p/1e6)
    return out
imp = {cas: {f"{T*1e6:.1f}": impact(T, J_eff(cas)) for T in T_RANGE[cas]} for cas in dp.CASSETTES}
beats_20h = F_BEAT * 20 * 3600
res["impact"] = dict(E_star_GPa=E_star/1e9, k_stop_N_per_m=K_STOP, per_case=imp, beats_per_20h=beats_20h,
                     petg_yield_MPa=50.0,
                     energy_ledger_uJ={f"{T*1e6:.1f}": dict(unlatch=float(np.interp(0.2, MU, E_unl))*1e6, impulse=E_imp*1e6,
                                                           cocking_losses=(TAU_ARM*DTH_ARM/dp.LIFT_EFF - TAU_ARM*DTH_ARM)*1e6,
                                                           lock_impact=impact(T, J_eff("10 mL/h"))["E_surplus_uJ"])
                                       for T in (19.5e-6, 100e-6, 476e-6)})
res["impact"]["worst"] = imp["10 mL/h"][f"{T_RANGE['10 mL/h'][1]*1e6:.1f}"]

# ------------------------------------------------------------------ 6. rate sensitivity via the beat-level simulation
I, k = 5.178802761191868e-06, 0.00029523676096236847   # from results/system_results.json (oscillator)
try:
    _r = json.load(open("results/system_results.json"))["oscillator"]; I, k = _r["I_kgm2"], _r["k_Nm_per_rad"]
except Exception:
    pass
base_E, base_U = E_imp, dp.E_UNLATCH                 # rev 5: base impulse = forward-model value
sens = []
for fE in (0.9, 1.0, 1.1):
    for U in (E_unl[0], E_unl[1], E_unl[3]):
        dp.E_IMPULSE, dp.E_UNLATCH = base_E * fE, float(U)
        r = ae.simulate(I, k, 150.0, "arm", t_end=120.0)
        sens.append(dict(E_impulse_factor=fE, E_unlatch_uJ=float(U)*1e6, status=r.get("status"),
                         rate_error_pct=100*r.get("rate_error", float("nan")), amplitude_rad=r.get("amplitude", float("nan")),
                         period_jitter_pct=100*r.get("period_jitter", float("nan"))))
dp.E_IMPULSE, dp.E_UNLATCH = base_E, base_U
res["rate_sensitivity"] = sens
res["rate_sensitivity_base_E_uJ"] = base_E*1e6
ok = [s for s in sens if s["status"] == "ok"]
res["rate_sensitivity_summary"] = dict(rate_error_min_pct=min(s["rate_error_pct"] for s in ok), rate_error_max_pct=max(s["rate_error_pct"] for s in ok),
                                       spread_pct=max(s["rate_error_pct"] for s in ok) - min(s["rate_error_pct"] for s in ok),
                                       amplitude_range_rad=[min(s["amplitude_rad"] for s in ok), max(s["amplitude_rad"] for s in ok)],
                                       n_ok=len(ok), n=len(sens))
json.dump(res, open("results/arm_mechanism.json", "w"), indent=1, default=float)

# ------------------------------------------------------------------ figure 7
apply_figure_style(sizes=(8, 7, 6))
fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.4))
ax = axs[0, 0]
cols = {"10 mL/h": "#2f6f9f", "5 mL/h": "#4f9f6f", "2 mL/h": "#e8a33d", "1 mL/h": "#c0392b"}
for cas, d in dyn.items():
    ax.plot(d["T_uNm"], d["t_recock_s"], color=cols[cas], lw=0.8, alpha=0.5)
    lo, hi = d["T_range_uNm"]; sel = (np.array(d["T_uNm"]) >= lo*0.98) & (np.array(d["T_uNm"]) <= hi*1.02)
    ax.plot(np.array(d["T_uNm"])[sel], np.array(d["t_recock_s"])[sel], color=cols[cas], lw=2.2, label=f"{cas} cassette, J = {d['J_eff_kgm2']*1e9:.0f} g·mm²")
ax.axhline(BUDGET, color=META_GREY, lw=0.8, ls="--"); ax.text(10.5, BUDGET*0.80, f"budget: half period minus impulse ({BUDGET:.2f} s)", ha="left", va="top", color=META_GREY, fontsize=6)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Escape-wheel torque (µN·m)"); ax.set_ylabel("Re-cock + lock time (s)")
ax.set_title("thick: each cassette's own torque range", fontsize=6, color=META_GREY, loc="left")
ax.legend(loc="lower left", fontsize=5.5); panel_letter(ax, "a")
ax = axs[0, 1]
Ts = np.array(dyn["10 mL/h"]["T_uNm"]) * 1e-6
for cas, ls in (("10 mL/h", "-"), ("1 mL/h", ":")):
    J = J_eff(cas)
    Es = [impact(T, J)["E_surplus_uJ"] for T in Ts]
    ax.plot(Ts*1e6, Es, color="k", ls=ls, lw=1.0, label=f"surplus energy per beat ({cas} train)")
ax.set_xscale("log"); ax.set_xlabel("Escape-wheel torque (µN·m)"); ax.set_ylabel("Energy into the lock per beat (µJ)")
ax2 = ax.twinx()
for name, k_s, c in ((list(K_STOP)[0], list(K_STOP.values())[0], "#c0392b"), (list(K_STOP)[1], list(K_STOP.values())[1], "#2f6f9f")):
    Fs = [math.sqrt(2 * max(impact(T, J_eff("10 mL/h"))["E_surplus_uJ"]*1e-6, 0) * k_s) for T in Ts]
    ax2.plot(Ts*1e6, Fs, color=c, lw=1.2, label=name.split(" (")[0])
ax2.set_yscale("log"); ax2.set_ylabel("Peak lock force (N)")
ax.axvspan(19.5, 476, color="#eeeeee", zorder=0)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=5.5, frameon=True, framealpha=0.92, edgecolor="none"); panel_letter(ax, "b")
ax = axs[1, 0]
labels = [f"{T:.0f}" for T in (19.5, 100, 476)]
led = res["impact"]["energy_ledger_uJ"]; keys = ["unlatch", "impulse", "cocking_losses", "lock_impact"]
kcol = {"unlatch": "#888888", "impulse": "#2f6f9f", "cocking_losses": "#e8a33d", "lock_impact": "#c0392b"}
bottom = np.zeros(3)
for key in keys:
    vals = np.array([led[k][key] for k in led]); ax.bar(labels, vals, bottom=bottom, color=kcol[key], label=key.replace("_", " "), width=0.6); bottom += vals
ax.set_xlabel("Escape-wheel torque (µN·m)"); ax.set_ylabel("Per-beat energy (µJ)"); ax.set_yscale("log"); ax.set_ylim(0.05, 300)
ax.legend(loc="upper left", fontsize=5.5); panel_letter(ax, "c")
ax = axs[1, 1]
for fE, m in ((0.9, "v"), (1.0, "o"), (1.1, "^")):
    rows = [s for s in sens if s["E_impulse_factor"] == fE and s["status"] == "ok"]
    ax.plot([s["E_unlatch_uJ"] for s in rows], [s["rate_error_pct"] for s in rows], marker=m, ms=4, lw=1.0, color="#2f6f9f",
            alpha=0.5 + 0.25*(fE-0.9)/0.2, label=f"E_impulse × {fE:.1f}")
ax.set_xlabel("Unlatch energy (µJ)  [µ = 0.1 … 0.4]"); ax.set_ylabel("Rate error vs free period (%)")
ax.legend(fontsize=5.5); panel_letter(ax, "d")
fig.tight_layout(); fig.savefig("figures/fig7_arm_mechanism.png", dpi=300); plt.close(fig)
print(json.dumps({k: res[k] for k in ("impulse", "unlatch", "lock")}, indent=0, default=float)[:1500])
print("recock:", {c: (round(d["T_range_uNm"][0],1), round(d["t_recock_at_min_torque_s"],3), round(d["time_margin_at_min_torque"],2), f"J={d['J_eff_kgm2']:.2e}") for c, d in dyn.items()})
print("impact:", {c: {T: (round(v["E_surplus_uJ"], 2), round(v[list(K_STOP)[0]]["peak_force_N"], 2), round(v[list(K_STOP)[1]]["peak_force_N"], 3), round(v[list(K_STOP)[0]]["hertz_stress_MPa"]), round(v[list(K_STOP)[1]]["hertz_stress_MPa"])) for T, v in d.items()} for c, d in imp.items()})
print("sens:", res["rate_sensitivity_summary"])
