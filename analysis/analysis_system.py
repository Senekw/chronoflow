"""System-level analysis of ChronoFlow: flow tables, escapement study, operating map,
accuracy vs flow-restriction pumps, Monte-Carlo error budget, structural & energy checks.
Writes results/*.json|csv and figures/*.png.  Figure helpers (apply_figure_style, panel_letter,
end_of_line_labels) come from the figure-style kernel plugin already loaded in this kernel.
"""
import os, json, importlib
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
import design_params as dp, analysis_oscillator as ao, analysis_escapement as ae
importlib.reload(dp); importlib.reload(ao); importlib.reload(ae)

os.makedirs("results", exist_ok=True); os.makedirs("figures", exist_ok=True)
R = {}

# ---------------------------------------------------------------- oscillator
k = ao.k_total(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0, dp.BAND_T0)
m_bal, I_bal = ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG,
                                  dp.N_TRIM, dp.M_TRIM, dp.R_TRIM, dp.HUB_M, dp.HUB_R)
# regulate: choose trim radius that hits F_TARGET
rt = np.linspace(dp.R_TRIM_MIN, dp.R_TRIM_MAX, 601)
fs = np.array([ao.frequency(k, ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG,
               dp.N_TRIM, dp.M_TRIM, r, dp.HUB_M, dp.HUB_R)[1]) for r in rt])
r_reg = rt[np.argmin(abs(fs - dp.F_TARGET))]
m_bal, I_bal = ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG, dp.N_TRIM, dp.M_TRIM, r_reg, dp.HUB_M, dp.HUB_R)
f0 = ao.frequency(k, I_bal)
orient = ao.orientation_study(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0, m_bal, I_bal)
R["oscillator"] = dict(k_Nm_per_rad=k, balance_mass_g=m_bal*1e3, I_kgm2=I_bal, f_Hz=f0, f_target=dp.F_TARGET,
    trim_radius_mm=r_reg*1e3, regulation_range_Hz=[fs.min(), fs.max()],
    regulation_range_pct=[(fs.min()/dp.F_TARGET-1)*100, (fs.max()/dp.F_TARGET-1)*100],
    orientation_Hz={kk: float(v) for kk, v in orient.items()},
    orientation_max_pct=max(abs(v/f0-1) for v in orient.values())*100,
    shear_stress_at_theta_max_MPa=ao.max_shear_stress(dp.BAND_T, dp.THETA_MAX, dp.BAND_LH)/1e6,
    fatigue_limit_MPa=ao.TAU_FATIGUE_STEEL/1e6,
    tempco_per_K=ao.temperature_coefficient(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0),
    band_sag_mm=ao.band_sag(m_bal, dp.BAND_LH, dp.BAND_T0)*1e3)

# ---------------------------------------------------------------- flow tables
rows = []
for cname, z in dp.CASSETTES.items():
    rc = dp.cassette_ratio(z)
    for sname, s in dp.SYRINGES.items():
        rows.append(dict(cassette=cname, ratio=rc, syringe=sname, flow_mL_h=dp.flow_mL_h(dp.F_TARGET, rc, s["area_mm2"]),
                         uL_per_beat=dp.flow_mL_h(dp.F_TARGET, rc, s["area_mm2"])*1000/(2*dp.F_TARGET*3600)))
flow_tab = pd.DataFrame(rows); flow_tab.to_csv("results/flow_table.csv", index=False)

# ---------------------------------------------------------------- escapement torque sweep
Te_grid = np.logspace(-6, np.log10(5e-4), 26)
sweep = []
for Te in Te_grid:
    for model in ("lever", "arm"):
        r = ae.simulate(I_bal, k, 150, model, T_e=Te, t_end=60)
        sweep.append(dict(T_e=Te, model=model, **{kk: (float(v) if isinstance(v, (float, np.floating)) else v) for kk, v in r.items()}))
sw = pd.DataFrame(sweep); sw.to_csv("results/escapement_sweep.csv", index=False)
# arm escapement vs Q
qrows = []
for Q in (40, 60, 100, 150, 250, 400, 800):
    r = ae.simulate(I_bal, k, Q, "arm", T_e=5e-5, t_end=80)
    qrows.append(dict(Q=Q, **{kk: float(v) if isinstance(v, (float, np.floating)) else v for kk, v in r.items()}))
qt = pd.DataFrame(qrows); qt.to_csv("results/arm_vs_Q.csv", index=False)
R["escapement"] = dict(T_min_arm_Nm=(dp.E_IMPULSE/dp.LIFT_EFF)/(np.pi/dp.Z_ESCAPE),
    arm_rate_error_pct_Q150=float(qt.loc[qt.Q == 150, "rate_error"].iloc[0]*100),
    arm_rate_error_range_pct=[float(qt.rate_error.min()*100), float(qt.rate_error.max()*100)],
    arm_amplitude_range_rad=[float(qt.amplitude.min()), float(qt.amplitude.max())])

# ---------------------------------------------------------------- operating map (stall margin): syringe x cassette x back-pressure
T_MIN = R["escapement"]["T_min_arm_Nm"]
P_GRID = (0.0, 100.0, 200.0, 300.0)                       # mmHg downstream pressure; 300 = occlusion trip for every syringe
op = []
for sname, s in dp.SYRINGES.items():
    for cname, z in dp.CASSETTES.items():
        Q20 = dp.flow_mL_h(dp.F_TARGET, dp.cassette_ratio(z))          # carriage speed is set by the cassette, not the syringe
        for Pmm in P_GRID:
            F_plunger = Pmm * 133.322 * s["area_mm2"] * 1e-6 + s["friction"]
            Fx = dp.F_SPRING - F_plunger
            Te, P, eta = ae.escape_wheel_torque(Fx, Q20)
            op.append(dict(syringe=sname, cassette=cname, pressure_mmHg=Pmm, F_plunger_N=F_plunger, F_excess_N=Fx,
                           power_uW=P*1e6, T_e_uNm=Te*1e6, margin=Te/T_MIN, flow_mL_h=dp.flow_mL_h(dp.F_TARGET, dp.cassette_ratio(z), s["area_mm2"])))
opm = pd.DataFrame(op); opm.to_csv("results/operating_map.csv", index=False)
F_min = float(opm.F_excess_N.min()); F_max = float(opm.F_excess_N.max())
worst = opm.loc[opm.margin.idxmin()]
# the same map for the rev-1 drive spring (30 N), to document why the spring was raised
op1 = opm.copy(); op1["F_excess_N"] = 30.0 - op1.F_plunger_N
op1["T_e_uNm"] = [ae.escape_wheel_torque(fx, dp.flow_mL_h(dp.F_TARGET, dp.cassette_ratio(dp.CASSETTES[c])))[0]*1e6 for fx, c in zip(op1.F_excess_N, op1.cassette)]
op1["margin"] = op1.T_e_uNm*1e-6 / T_MIN; op1.to_csv("results/operating_map_rev1_30N.csv", index=False)
R["operating"] = dict(F_spring_N=dp.F_SPRING, F_excess_min_N=F_min, F_excess_max_N=F_max, worst_margin=float(worst.margin),
                      worst_case=dict(syringe=worst.syringe, cassette=worst.cassette, pressure_mmHg=float(worst.pressure_mmHg)),
                      train_efficiency=float(eta), n_gear_meshes=len(dp.STAGES) + 2,
                      power_range_uW=[float(opm.power_uW.min()), float(opm.power_uW.max())],
                      torque_range_uNm=[float(opm.T_e_uNm.min()), float(opm.T_e_uNm.max())],
                      torque_dynamic_range=float(opm.T_e_uNm.max()/opm.T_e_uNm.min()),
                      rev1_30N=dict(worst_margin=float(op1.margin.min()), n_stall_cells=int((op1.margin < 1).sum()),
                                    n_below_1p7=int((op1.margin < 1.7).sum()), n_cells=int(len(op1)),
                                    stalled=[f"{r.syringe} / {r.cassette} / {r.pressure_mmHg:.0f} mmHg" for r in op1[op1.margin < 1].itertuples()]))

# ---------------------------------------------------------------- accuracy vs flow-restriction pumps
def mu_water(T_C):  # Vogel equation, mPa s
    T = T_C + 273.15
    return 0.02939 * np.exp(507.88 / (T - 149.3))
T_amb = np.linspace(10, 40, 61)
restr_T_err = (mu_water(20) / mu_water(T_amb) - 1) * 100            # calibrated at 20 C
chrono_T_err = R["oscillator"]["tempco_per_K"] * (T_amb - 20) * 100
Pback = np.linspace(0, 300, 61)
restr_P_err = {Pd: -Pback / Pd * 100 for Pd in (300.0, 500.0)}
chrono_P_err = np.where(Pback < 300, 0.0, np.nan)
R["comparison"] = dict(restriction_err_at_37C_pct=float((mu_water(20)/mu_water(37)-1)*100),
                       restriction_err_at_10C_pct=float((mu_water(20)/mu_water(10)-1)*100),
                       chrono_err_10_to_40C_pct=[float(chrono_T_err.min()), float(chrono_T_err.max())],
                       restriction_err_100mmHg_back_Pd300=float(-100/300*100))

# ---------------------------------------------------------------- Monte-Carlo steady-state flow error budget
rng = np.random.default_rng(7)
N = 20000
d_cal = rng.uniform(-0.3, 0.3, N)                                   # residual after trim regulation, %
Tamb = rng.uniform(15, 35, N); d_T = R["oscillator"]["tempco_per_K"] * (Tamb - 25) * 100
d_or = rng.uniform(-R["oscillator"]["orientation_max_pct"], R["oscillator"]["orientation_max_pct"], N)
Qs = rng.uniform(60, 400, N)
esc = np.interp(np.log(Qs), np.log(qt.Q), qt.rate_error*100)        # escapement error vs Q (rate error of PERIOD)
d_esc = -(esc - esc.mean())                                          # calibration removes the mean
d_syr = rng.uniform(-1.0, 1.0, N)                                    # syringe bore area tolerance (ISO 7886-1 class), %
d_rack = rng.uniform(-0.2, 0.2, N)                                   # residual rack pitch after print-scale calibration, %
total = d_cal + d_T + d_or + d_esc + d_syr + d_rack
mc = pd.DataFrame(dict(calibration=d_cal, temperature=d_T, orientation=d_or, escapement=d_esc, syringe=d_syr, rack=d_rack, total=total))
mc.describe().to_csv("results/montecarlo_summary.csv")
R["montecarlo"] = dict(p2_5=float(np.percentile(total, 2.5)), p97_5=float(np.percentile(total, 97.5)),
                       p99=float(np.percentile(abs(total), 99)), sd=float(total.std()),
                       contrib_sd_pct={c: float(mc[c].std()) for c in mc.columns if c != "total"})

# ---------------------------------------------------------------- structural / power / resolution / occlusion
Ft = dp.F_SPRING                      # tangential load on rack pinion teeth (worst: carriage fully loaded, static)
Y12 = 0.245                           # Lewis form factor, 12 T, 20 deg
sigma_stage1 = Ft / (dp.GEAR_FACE_STAGE1 * dp.MODULE * Y12)         # MPa
torques = []
T_shaft = dp.F_SPRING * (dp.MODULE * dp.Z_RACK_PINION / 2) * 1e-3   # N m on rack pinion shaft
for i, (w, p) in enumerate(dp.STAGES):
    torques.append(dict(stage=i+1, torque_in_mNm=T_shaft*1e3))
    T_shaft = T_shaft / (w / p) * dp.GEAR_EFF
R["structure"] = dict(stage1_tooth_bending_MPa=sigma_stage1, petg_yield_MPa=50.0, safety_factor=50.0/sigma_stage1,
                      rack_pinion_torque_mNm=torques[0]["torque_in_mNm"], stage_torques=torques)
E_spring = dp.F_SPRING * dp.TRAVEL * 1e-3
E_esc_per_beat = dp.E_IMPULSE / dp.LIFT_EFF + dp.E_UNLATCH
P_esc = E_esc_per_beat * 2 * dp.F_TARGET
R["energy"] = dict(spring_energy_J=E_spring, escapement_power_uW=P_esc*1e6,
                   escapement_energy_20h_J=P_esc*20*3600,
                   runtime_h_20mL={c: 20/dp.flow_mL_h(dp.F_TARGET, dp.cassette_ratio(z)) for c, z in dp.CASSETTES.items()})
occ = {}
for sname, s in dp.SYRINGES.items():
    F_set = dp.trip_force(sname); F_pre = dp.trip_preload(sname)
    lo, nom, hi = dp.trip_pressure_range(sname)
    A = s["area_mm2"]
    occ[sname] = dict(trip_force_N=F_set, spacer_preload_N=F_pre, spacer_preload_mm=F_pre/dp.TRIP_SPRING_RATE,
                      sense_travel_force_rise_N=dp.TRIP_SPRING_RATE*dp.TRIP_SENSE_TRAVEL,
                      sense_rise_mmHg_if_uncompensated=dp.TRIP_SPRING_RATE*dp.TRIP_SENSE_TRAVEL/(A*1e-6)/133.322,
                      full_3mm_rise_mmHg=dp.TRIP_SPRING_RATE*3.0/(A*1e-6)/133.322,
                      pressure_band_mmHg=[lo, nom, hi], pressure_band_pct=[(lo/nom-1)*100, (hi/nom-1)*100],
                      stored_displacement_at_sear_mL=dp.stored_displacement_mL(sname), stored_displacement_3mm_mL=A*3.0/1000.0)
R["occlusion"] = dict(per_syringe=occ, spring_rate_N_per_mm=dp.TRIP_SPRING_RATE, sense_travel_mm=dp.TRIP_SENSE_TRAVEL, flag_travel_mm=dp.FLAG_TRAVEL,
                      preload_range_N=[min(v["spacer_preload_N"] for v in occ.values()), max(v["spacer_preload_N"] for v in occ.values())],
                      friction_uncertainty_N=dp.FRICTION_UNCERTAINTY)
# band adjustment: symmetric vs one-sided coverage of ribbon thickness tolerance
from scipy.optimize import brentq as _brentq
def _f_lr(Lu, Ll, t): return ao.frequency(ao.k_half(dp.BAND_W, t, Lu, dp.BAND_T0) + ao.k_half(dp.BAND_W, t, Ll, dp.BAND_T0), I_bal)
_f0 = _f_lr(dp.BAND_LH, dp.BAND_LH, dp.BAND_T)
band_adj = {}
for tag, t in (("t_min", dp.BAND_T - dp.BAND_T_TOL), ("t_max", dp.BAND_T + dp.BAND_T_TOL)):
    band_adj[tag] = dict(t_um=t*1e6, L_symmetric_mm=_brentq(lambda L: _f_lr(L, L, t) - _f0, 5e-3, 80e-3)*1e3,
                         L_one_sided_mm=_brentq(lambda L: _f_lr(L, dp.BAND_LH, t) - _f0, 5e-3, 80e-3)*1e3)
band_adj["one_sided_coverage_um"] = [_brentq(lambda t: _f_lr(L, dp.BAND_LH, t) - _f0, 20e-6, 60e-6)*1e6 for L in (dp.BAND_LH_MIN, dp.BAND_LH_MAX)]
band_adj["sensitivity_pct_per_mm"] = dict(symmetric=(_f_lr(dp.BAND_LH+1e-3, dp.BAND_LH+1e-3, dp.BAND_T)/_f0-1)*100, one_sided=(_f_lr(dp.BAND_LH+1e-3, dp.BAND_LH, dp.BAND_T)/_f0-1)*100)
R["oscillator"]["band_adjustment"] = band_adj
R["resolution_uL_per_beat"] = {c: float(flow_tab[(flow_tab.cassette == c) & (flow_tab.syringe == dp.DESIGN_SYRINGE)].uL_per_beat.iloc[0]) for c in dp.CASSETTES}

# ---------------------------------------------------------------- BOM (indicative unit prices, USD, small quantities)
bom = [
 ("Printed parts (PETG, ~420 g incl. supports)", 1, 9.00),
 ("Constant-force spring 5 lbf (22.2 N) x 1.5 m (x2, parallel)", 2, 7.50),
 ("Feeler gauge 0.0015\" (38.1 um) blade (taut band stock)", 1, 1.50),
 ("Steel coil spring, tensioner (1 N @ 5 mm)", 1, 0.40),
 ("Steel pins 3 mm x 50 mm (gear shafts)", 8, 0.30),
 ("Steel pins 1 mm (impulse pin, pallets, latches)", 4, 0.10),
 ("Brass M3 hex nuts (trim weights)", 4, 0.08),
 ("Compression spring, occlusion trip 3 N/mm + 5 keyed preload spacers (5-26 N)", 1, 0.90),
 ("Flag-rod spring 2 N + steel sear pin (snap-action trip)", 1, 0.40),
 ("M3 hardware (screws, heat-set inserts)", 20, 0.08),
 ("Linear rods 6 mm x 200 mm (carriage rails)", 2, 1.80),
 ("Rate cassettes (printed gear sets), 4", 4, 1.20),
]
bom_df = pd.DataFrame(bom, columns=["item", "qty", "unit_usd"]); bom_df["ext_usd"] = bom_df.qty*bom_df.unit_usd
bom_df.to_csv("results/bom.csv", index=False)
R["bom_total_usd"] = float(bom_df.ext_usd.sum())

json.dump(R, open("results/system_results.json", "w"), indent=1, default=float)

# ================================================================ FIGURES
apply_figure_style(sizes=(8, 7, 6))
C_ARM, C_LEV, C_RESTR = "#1b6ca8", "#b05a1a", "#7a7a7a"

# --- Fig 2: escapement comparison
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
ax = axs[0]
for model, col, lab in (("lever", C_LEV, "conventional lever escapement"), ("arm", C_ARM, "constant-force spring-arm (ChronoFlow)")):
    d = sw[sw.model == model]
    ok = d[d.status == "ok"]
    ax.plot(ok.T_e*1e6, ok.amplitude, "-o", ms=3, color=col, label=lab)
    ob = d[d.status == "overbank"]
    if len(ob): ax.plot(ob.T_e*1e6, np.full(len(ob), dp.THETA_MAX), "x", color=col, ms=5)
    st = d[d.status == "stall"]
    if len(st):
        ax.plot(st.T_e*1e6, np.zeros(len(st)), "|", color=col, ms=8)
        if model == "arm": ax.text(st.T_e.iloc[0]*1e6, 0.1, "arm not re-cocked\n(below minimum torque)", fontsize=6, color=col, va="bottom")
ax.axhline(dp.THETA_MAX, color=META_GREY, lw=0.8, ls="--"); ax.text(480, dp.THETA_MAX+0.03, "banking limit", color=META_GREY, fontsize=7, ha="right")
lo, hi = opm.T_e_uNm.min(), opm.T_e_uNm.max()
ax.axvspan(lo, hi, color="#dddddd", alpha=0.5, lw=0); ax.text(np.sqrt(lo*hi), 0.08, "operating range\n(all rates, all loads)", ha="center", fontsize=7, color="#444444")
ax.set_xscale("log"); ax.set_xlabel("escape-wheel torque (µN·m)"); ax.set_ylabel("balance amplitude (rad)")
ax.set_title(f"Spring-arm: constant amplitude over {R['operating']['torque_dynamic_range']:.0f}× torque")
ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"]); ax.set_ylim(-0.05, 1.8); ax.margins(x=0.04)
ax.legend(frameon=False, loc="upper left", fontsize=7, bbox_to_anchor=(0.0, 0.93))
ax = axs[1]
for model, col in (("lever", C_LEV), ("arm", C_ARM)):
    ok = sw[(sw.model == model) & (sw.status == "ok")]
    ax.plot(ok.T_e*1e6, ok.rate_error*100, "-o", ms=3, color=col)
ax.axvspan(lo, hi, color="#dddddd", alpha=0.5, lw=0)
ax.axhline(0, color=META_GREY, lw=0.6)
ax.set_xscale("log"); ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
ax.set_xlabel("escape-wheel torque (µN·m)"); ax.set_ylabel("period error vs free oscillator (%)")
ax.set_title("Escapement rate error < 0.1 % in operating range"); ax.margins(0.08)
panel_letter(axs[0], "a"); panel_letter(axs[1], "b")
fig.savefig("figures/fig2_escapement.png", dpi=300)

# --- Fig 3: accuracy vs restriction pumps
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)
ax = axs[0]
ax.plot(T_amb, restr_T_err, color=C_RESTR, lw=1.6)
ax.plot(T_amb, chrono_T_err, color=C_ARM, lw=2.0)
ax.axhspan(-5, 5, color="#e8f1e8", lw=0); ax.text(10.5, 4.0, "±5 % clinical tolerance", fontsize=7, color="#3a6b3a")
ax.axhline(0, color=META_GREY, lw=0.6)
ax.set_xlabel("fluid temperature (°C)"); ax.set_ylabel("flow-rate error (%)")
ax.set_title("Volumetric metering removes the viscosity dependence")
ax.text(T_amb[-1]+0.5, restr_T_err[-1], "flow-restriction pump\n(Q ∝ ΔP/μ)", color=C_RESTR, fontsize=7, va="center")
ax.text(T_amb[-1]+0.5, chrono_T_err[-1]-3, "ChronoFlow\n(volumetric)", color=C_ARM, fontsize=7, va="center")
ax.margins(x=0.04); ax.set_xlim(10, 47)
ax = axs[1]
ax.plot(Pback, restr_P_err[300.0], color=C_RESTR, lw=1.6, label="restriction pump, 300 mmHg drive")
ax.plot(Pback, restr_P_err[500.0], color=C_RESTR, lw=1.2, ls="--", label="restriction pump, 500 mmHg drive")
ax.plot(Pback, chrono_P_err, color=C_ARM, lw=2.0, label="ChronoFlow (stops + flags at 300 mmHg)")
ax.plot([300], [0], "v", color="#c0392b", ms=6); ax.text(300, 4, "occlusion trip", ha="right", fontsize=7, color="#c0392b")
ax.axhspan(-5, 5, color="#e8f1e8", lw=0); ax.axhline(0, color=META_GREY, lw=0.6)
ax.set_xlabel("downstream back-pressure (mmHg)"); ax.set_ylabel("flow-rate error (%)")
ax.set_title("Back-pressure does not change the metered rate"); ax.legend(frameon=False, fontsize=7, loc="lower left"); ax.margins(0.05)
panel_letter(axs[0], "a"); panel_letter(axs[1], "b")
fig.savefig("figures/fig3_accuracy_comparison.png", dpi=300)

# --- Fig 4: error budget
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.1, 1]))
ax = axs[0]
contrib = pd.Series(R["montecarlo"]["contrib_sd_pct"]).sort_values()
names = {"syringe": "syringe bore tolerance (±1 %)", "calibration": "trim regulation residual (±0.3 %)", "rack": "rack pitch after print calibration (±0.2 %)",
         "temperature": "ambient 15–35 °C", "orientation": "device orientation", "escapement": "escapement error, Q = 60–400"}
ax.barh([names[c] for c in contrib.index], contrib.values, color=C_ARM, height=0.6)
for y, v in enumerate(contrib.values): ax.text(v + 0.01, y, f"{v:.2f}", va="center", fontsize=6)
ax.set_xlabel("standard deviation of flow error (%)"); ax.set_title("Syringe bore tolerance dominates"); ax.margins(x=0.15)
ax = axs[1]
ax.hist(total, bins=80, color=C_ARM, alpha=0.85)
p = R["montecarlo"]
ax.axvline(p["p2_5"], color=META_GREY, ls="--", lw=0.8); ax.axvline(p["p97_5"], color=META_GREY, ls="--", lw=0.8)
ax.text(0.02, 0.95, f"95 % of units within\n{p['p2_5']:+.2f} … {p['p97_5']:+.2f} %\n(n = {N:,} Monte-Carlo draws)", transform=ax.transAxes, va="top", fontsize=7)
ax.set_xlabel("steady-state flow-rate error (%)"); ax.set_ylabel("count"); ax.set_title("Predicted accuracy is well inside ±5 %")
ax.set_xlim(-5, 5); ax.axvspan(-5, 5, color="#e8f1e8", lw=0, zorder=0)
panel_letter(axs[0], "a"); panel_letter(axs[1], "b")
fig.savefig("figures/fig4_error_budget.png", dpi=300)

# --- Fig 5: operating map (stall margin over syringe x cassette, at the 300 mmHg trip and free-running)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.15, 1]))
cas = list(dp.CASSETTES); syr = list(dp.SYRINGES)
ax = axs[0]
M = np.array([[opm[(opm.syringe == s_) & (opm.cassette == c_) & (opm.pressure_mmHg == 300.0)].margin.iloc[0] for s_ in syr] for c_ in cas])
im = ax.imshow(M, cmap="Blues", norm=mpl.colors.LogNorm(vmin=1, vmax=M.max()*1.1), aspect="auto")
for i in range(len(cas)):
    for j in range(len(syr)):
        ax.text(j, i, f"{M[i, j]:.1f}×", ha="center", va="center", fontsize=7, color="white" if M[i, j] > 4 else "black")
ax.set_xticks(range(len(syr))); ax.set_xticklabels([s_.replace("BD ", "") for s_ in syr]); ax.set_yticks(range(len(cas))); ax.set_yticklabels([c_.replace(" mL/h", "") for c_ in cas])
ax.set_xlabel("syringe (BD)"); ax.set_ylabel("cassette (mL/h on 20 mL)")
ax.set_title(f"Margin at the 300 mmHg trip (worst {M.min():.1f}×: {worst.syringe.replace('BD ', '')}, {worst.cassette})")
ax = axs[1]
x = np.arange(len(cas))
for j, (Pmm, col, lab) in enumerate(((300.0, "#9ecae1", "60 mL at 300 mmHg trip (excess %.1f N)" % F_min), (0.0, C_ARM, "5 mL free-running (excess %.1f N)" % F_max))):
    sname = "BD 60 mL" if Pmm == 300.0 else "BD 5 mL"
    d = opm[(opm.syringe == sname) & (opm.pressure_mmHg == Pmm)].set_index("cassette").loc[cas]
    ax.bar(x + (j-0.5)*0.38, d.margin, width=0.36, color=col, label=lab)
ax.axhline(1, color="#c0392b", lw=0.9); ax.text(len(cas)-0.55, 1.12, "stall threshold", color="#c0392b", fontsize=7, ha="right")
r1 = op1[(op1.syringe == "BD 60 mL") & (op1.pressure_mmHg == 300.0)].set_index("cassette").loc[cas]
ax.plot(x - 0.19, r1.margin, "x", color="#c0392b", ms=6, label="rev 1 (30 N spring): 60 mL at trip")
ax.set_yscale("log"); ax.set_yticks([0.3, 1, 10, 100]); ax.set_yticklabels(["0.3×", "1×", "10×", "100×"])
ax.set_xticks(x); ax.set_xticklabels([c_.replace(" mL/h", "") for c_ in cas]); ax.set_xlabel("rate cassette (mL/h, BD 20 mL)")
ax.set_ylabel("escape-wheel torque / re-cock minimum"); ax.set_title(f"Extremes of the map: {F_min:.0f}–{F_max:.0f} N excess")
ax.legend(frameon=False, fontsize=6.5, loc="upper right"); ax.margins(y=0.15)
panel_letter(axs[0], "a"); panel_letter(axs[1], "b")
fig.savefig("figures/fig5_operating_map.png", dpi=300)
print(json.dumps({kk: R[kk] for kk in ("oscillator", "escapement", "operating", "montecarlo", "comparison", "structure", "energy", "bom_total_usd", "resolution_uL_per_beat")}, indent=0, default=float)[:3500])
