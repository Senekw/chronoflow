"""Render the data/schematic panels of the manuscript figures at exact slot pixel sizes, then compose.
Requires the figure-style and figure-composer kernel plugins (apply_figure_style, panel_letter,
compose_figure, META_GREY) to be loaded in the kernel that exec()s this file."""
import json, importlib, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Polygon, FancyBboxPatch, Wedge
import design_params as dp, analysis_oscillator as ao, analysis_escapement as ae, analysis_crossspring as acs
for m in (dp, ao, ae, acs): importlib.reload(m)

DPI = 300
from matplotlib import font_manager as _fm
FONT_BOLD = _fm.findfont(_fm.FontProperties(family="DejaVu Sans", weight="bold"))
OUT = json.load(open("panels/outlines.json")); PX = json.load(open("panels/panel_px.json"))
apply_figure_style(sizes=(7.5, 6.5, 6))
C_ARM, C_LEV, C_RESTR, C_ALARM, C_OK = "#1b6ca8", "#b05a1a", "#7a7a7a", "#c0392b", "#e8f1e8"
C_STEEL, C_PRINT, C_BRASS = "#8c96a0", "#2f6f9f", "#c8a951"

def new(fig_key, L, **kw):
    w, h = PX[fig_key][L]
    fig, ax = plt.subplots(figsize=(w/DPI, h/DPI), dpi=DPI, **kw)
    return fig, ax
def save(fig, fig_key, L):
    fig.savefig(f"panels/{fig_key}_{L}.png", dpi=DPI, facecolor="white")
    plt.close(fig)

# ---------------------------------------------------------------- shared physics
k = ao.k_total(dp.BAND_W, dp.BAND_T, dp.BAND_LH, dp.BAND_T0, dp.BAND_T0)
R = json.load(open("results/system_results.json"))
r_reg = R["oscillator"]["trim_radius_mm"] * 1e-3
m_bal, I_bal = ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG, dp.N_TRIM, dp.M_TRIM, r_reg, dp.HUB_M, dp.HUB_R)
sw = pd.read_csv("results/escapement_sweep.csv"); qt = pd.read_csv("results/arm_vs_Q.csv"); opm = pd.read_csv("results/operating_map.csv")
ft = pd.read_csv("results/flow_table.csv")
lo_T, hi_T = opm.T_e_uNm.min(), opm.T_e_uNm.max()

def mu_water(T_C):
    T = T_C + 273.15; return 0.02939 * np.exp(507.88 / (T - 149.3))

# ================================================================ FIG 1
# b: block diagram of the energy / metering chain
fig, ax = new("fig1", "b"); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10.6)
boxes = [(f"constant-force\nspring, {dp.F_SPRING:.1f} N", 9.6, C_STEEL), ("carriage + pusher\n(plunger load 1.5–27 N)", 8.1, C_PRINT), ("rack → pinion\n+ winding ratchet", 6.6, C_PRINT),
         ("step-up train\n5·5·5·2.5 = 312.5", 5.1, C_PRINT), ("rate cassette\n×1 · ×2 · ×5 · ×10", 3.6, "#e8a33d"), ("escape wheel 15 T +\nspring-arm escapement", 2.1, C_BRASS),
         ("taut-band balance\nf₀ = 1.20 Hz", 0.6, C_PRINT)]
for txt, y, col in boxes:
    ax.add_patch(FancyBboxPatch((1.2, y-0.5), 5.6, 1.0, boxstyle="round,pad=0.05", fc="white", ec=col, lw=1.4))
    ax.text(4.0, y, txt, ha="center", va="center", fontsize=6.5)
for y in [b[1] for b in boxes[:-1]]:
    ax.add_patch(FancyArrowPatch((4.0, y-0.5), (4.0, y-1.0), arrowstyle="-|>", mutation_scale=7, lw=1.0, color="#333333"))
# side annotations
ax.text(7.2, 9.6, f"F = const\n({R['energy']['spring_energy_J']:.1f} J per wind)", fontsize=6, va="center", color="#444444")
ax.text(7.2, 8.1, f"excess force\n{R['operating']['F_excess_min_N']:.0f}–{R['operating']['F_excess_max_N']:.0f} N", fontsize=6, va="center", color="#444444")
ax.text(7.2, 6.6, f"{R['structure']['rack_pinion_torque_mNm']:.0f} mN·m", fontsize=6, va="center", color="#444444")
ax.text(7.2, 5.1, "torque ÷ 312", fontsize=6, va="center", color="#444444")
ax.text(7.2, 3.6, "sets rate", fontsize=6, va="center", color="#444444")
ax.text(7.2, 2.1, f"{R['operating']['torque_range_uNm'][0]:.0f}–{R['operating']['torque_range_uNm'][1]:.0f} µN·m in,\nfixed 1.2 µJ out", fontsize=6, va="center", color="#444444")
ax.text(7.2, 0.6, "1 half-tooth / beat\n→ 4 µm, 1.16 µL", fontsize=6, va="center", color="#444444")
ax.text(0.15, 5.1, "energy", rotation=90, va="center", fontsize=6.5, color=C_STEEL)
ax.add_patch(FancyArrowPatch((0.6, 9.3), (0.6, 2.6), arrowstyle="-|>", mutation_scale=7, lw=1.0, color=C_STEEL))
ax.text(9.7, 5.1, "rate", rotation=270, va="center", fontsize=6.5, color=C_ARM)
ax.add_patch(FancyArrowPatch((9.4, 1.0), (9.4, 7.8), arrowstyle="-|>", mutation_scale=7, lw=1.0, color=C_ARM))
save(fig, "fig1", "b")

# c: temperature
T_amb = np.linspace(10, 40, 61)
fig, ax = new("fig1", "c")
ax.axhspan(-5, 5, color=C_OK, lw=0); ax.text(10.5, 6.5, "±5 % clinical tolerance", fontsize=6, color="#3a6b3a")
ax.plot(T_amb, (mu_water(20)/mu_water(T_amb)-1)*100, color=C_RESTR, lw=1.5)
ax.plot(T_amb, R["oscillator"]["tempco_per_K"]*(T_amb-20)*100, color=C_ARM, lw=2.0)
ax.axhline(0, color=META_GREY, lw=0.6)
ax.text(40.6, 44, "flow-restriction\npump (Q ∝ ΔP/μ)", color=C_RESTR, fontsize=6.5, va="center")
ax.text(40.6, -1, "ChronoFlow\n(volumetric)", color=C_ARM, fontsize=6.5, va="center")
ax.set_xlabel("fluid temperature (°C)"); ax.set_ylabel("flow-rate error (%)"); ax.set_xlim(10, 49); ax.set_ylim(-28, 58)
ax.set_title("Viscosity no longer sets the rate", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig1", "c")

# d: back-pressure
Pb = np.linspace(0, 300, 61)
fig, ax = new("fig1", "d")
ax.axhspan(-5, 5, color=C_OK, lw=0); ax.axhline(0, color=META_GREY, lw=0.6)
ax.plot(Pb, -Pb/300*100, color=C_RESTR, lw=1.5); ax.plot(Pb, -Pb/500*100, color=C_RESTR, lw=1.1, ls="--")
ax.plot(Pb, np.zeros_like(Pb), color=C_ARM, lw=2.0)
ax.plot([300], [0], "v", color=C_ALARM, ms=5); ax.text(306, 2, "occlusion trip\n(stops + flags)", ha="left", fontsize=6, color=C_ALARM, va="bottom")
ax.text(305, -100, "restriction pump,\n300 mmHg drive", color=C_RESTR, fontsize=6.5, va="center")
ax.text(305, -60, "restriction pump,\n500 mmHg drive", color=C_RESTR, fontsize=6.5, va="center")
ax.text(120, 7, "ChronoFlow (volumetric)", color=C_ARM, fontsize=6.5, va="bottom")
ax.set_xlabel("downstream back-pressure (mmHg)"); ax.set_ylabel("flow-rate error (%)"); ax.set_xlim(0, 380); ax.set_ylim(-110, 22)
ax.set_title("Back-pressure no longer sets the rate", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig1", "d")

# ================================================================ FIG 2
# a: escapement mechanism schematic (three phases)
fig, ax = new("fig2", "a"); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 7.6); ax.set_aspect("equal"); fig.subplots_adjust(0, 0, 1, 1)
def draw_phase(x0, phase):
    # escape wheel
    cx, cy = x0 + 1.3, 4.9
    th = np.linspace(0, 2*np.pi, 16)[:-1]
    pts = []
    for a in th:
        pts += [(cx + 1.1*np.cos(a), cy + 1.1*np.sin(a)), (cx + 0.8*np.cos(a+0.12), cy + 0.8*np.sin(a+0.12)), (cx + 0.8*np.cos(a+0.34), cy + 0.8*np.sin(a+0.34))]
    ax.add_patch(Polygon(pts, closed=True, fc=C_BRASS, ec="#7a6230", lw=0.6))
    ax.add_patch(Circle((cx, cy), 0.12, fc="white", ec="#333"))
    # arm pivot and arm
    px, py = x0 + 1.5, 3.2
    ang = {"cock": 0.0, "latch": 0.0, "release": -0.35}[phase]
    L = 2.0
    ex, ey = px + L*np.cos(0.35+ang), py + L*np.sin(0.35+ang)
    ax.plot([px, ex], [py, ey], color="#e8a33d", lw=3, solid_capstyle="round")
    ax.add_patch(Circle((px, py), 0.12, fc="white", ec="#333", zorder=3))
    # pallet pin on arm toward wheel
    ax.plot([px, px+0.15], [py, py+1.05], color="#e8a33d", lw=3, solid_capstyle="round")
    ax.plot(px+0.15, py+1.05, "o", ms=3, color=C_STEEL)
    # flexure spring (zigzag) from frame post to arm
    zx = np.linspace(px-0.9, px-0.1, 9); zy = py - 0.7 + 0.12*np.array([1,-1]*4+[1])
    ax.plot(zx, zy, color="#555", lw=0.8); ax.add_patch(Rectangle((px-1.05, py-0.95), 0.15, 0.5, fc=C_STEEL, ec="none"))
    # latch
    lx, ly = ex - 0.1, ey + 0.25
    if phase != "release":
        ax.add_patch(Rectangle((lx-0.1, ly), 0.22, 0.5, fc=C_STEEL, ec="#333", lw=0.5))
    else:
        ax.add_patch(Rectangle((lx+0.2, ly+0.3), 0.22, 0.5, fc=C_STEEL, ec="#333", lw=0.5, alpha=0.6))
    # balance hub + impulse pin
    bx, by = x0 + 3.5, 2.2
    ax.add_patch(Circle((bx, by), 0.55, fc="none", ec=C_PRINT, lw=1.8))
    pa = {"cock": 2.6, "latch": 2.0, "release": 1.5}[phase]
    ipx, ipy = bx + 0.75*np.cos(pa), by + 0.75*np.sin(pa)
    ax.plot(ipx, ipy, "o", ms=3.2, color=C_STEEL, zorder=4)
    ax.add_patch(FancyArrowPatch((bx+0.95*np.cos(pa+0.5), by+0.95*np.sin(pa+0.5)), (bx+0.95*np.cos(pa-0.5), by+0.95*np.sin(pa-0.5)),
                                 connectionstyle="arc3,rad=0.3", arrowstyle="-|>", mutation_scale=6, lw=0.8, color=C_PRINT))
    if phase == "release":
        ax.add_patch(FancyArrowPatch((ex-0.3, ey-0.2), (ex+0.25, ey-0.55), arrowstyle="-|>", mutation_scale=7, lw=1.2, color=C_ALARM))
    if phase == "cock":
        ax.add_patch(FancyArrowPatch((cx+0.3, cy-1.4), (cx-0.5, cy-1.3), connectionstyle="arc3,rad=0.4", arrowstyle="-|>", mutation_scale=6, lw=0.9, color="#7a6230"))
for i, (ph, t1, t2) in enumerate((("cock", "1  cock", "tooth lifts arm,\nflexure stores E"), ("latch", "2  hold", "latch holds arm;\nexcess torque idle"), ("release", "3  release", "pin trips latch;\nfixed 1.2 µJ impulse"))):
    draw_phase(4.0*i, ph)
    ax.text(4.0*i + 0.7, 7.3, t1, fontsize=7, va="top", fontweight="bold"); ax.text(4.0*i + 0.1, 6.75, t2, fontsize=5.8, va="top", color="#333")
ax.text(0.1, 0.3, "escape wheel", fontsize=5.5, color="#7a6230"); ax.text(2.6, 0.3, "arm · flexure · latch", fontsize=5.5, color="#9a6a1a"); ax.text(6.4, 0.3, "balance hub + impulse pin", fontsize=5.5, color=C_PRINT)
save(fig, "fig2", "a")

# b: time traces
fig, axs = plt.subplots(2, 1, figsize=(PX["fig2"]["b"][0]/DPI, PX["fig2"]["b"][1]/DPI), dpi=DPI, sharex=True)
for ax, model, col, lab in ((axs[0], "lever", C_LEV, "conventional lever escapement"), (axs[1], "arm", C_ARM, "constant-force spring-arm (ChronoFlow)")):
    for Te, ls, tl in ((lo_T*1e-6, "-", f"{lo_T:.0f} µN·m"), (hi_T*1e-6, ":", f"{hi_T:.0f} µN·m")):
        r = ae.simulate(I_bal, k, 150, model, T_e=Te, t_end=14, theta0=0.5, trace=True)
        ax.plot(r["t"], r["theta"], color=col, lw=0.9, ls=ls, label=f"T_e = {tl}" + ("  → over-banks" if r["status"] == "overbank" else ""))
    ax.axhline(dp.THETA_MAX, color=META_GREY, lw=0.6, ls="--"); ax.axhline(-dp.THETA_MAX, color=META_GREY, lw=0.6, ls="--")
    ax.set_ylim(-2.0, 2.0); ax.set_ylabel("θ (rad)")
    ax.text(0.01, 0.97, lab, transform=ax.transAxes, va="top", fontsize=6.5, color=col)
    ax.legend(frameon=False, fontsize=5.5, loc="lower left", ncol=2, bbox_to_anchor=(0.0, -0.02))
axs[0].text(13.8, dp.THETA_MAX+0.08, "banking limit", fontsize=5.5, color=META_GREY, ha="right")
axs[1].set_xlabel("time (s)"); axs[0].set_title("Lever amplitude follows torque; spring-arm amplitude does not", loc="left")
fig.tight_layout(pad=0.3, h_pad=0.4); save(fig, "fig2", "b")

# c: amplitude vs torque
fig, ax = new("fig2", "c")
for model, col in (("lever", C_LEV), ("arm", C_ARM)):
    d = sw[sw.model == model]; ok = d[d.status == "ok"]
    ax.plot(ok.T_e*1e6, ok.amplitude, "-o", ms=2.5, color=col)
    ob = d[d.status == "overbank"]
    if len(ob): ax.plot(ob.T_e*1e6, np.full(len(ob), dp.THETA_MAX), "x", color=col, ms=4)
    st = d[d.status == "stall"]
    if len(st): ax.plot(st.T_e*1e6, np.zeros(len(st)), "|", color=col, ms=6)
ax.axvspan(lo_T, hi_T, color="#e4e4e4", lw=0); ax.axhline(dp.THETA_MAX, color=META_GREY, lw=0.6, ls="--")
ax.text(np.sqrt(lo_T*hi_T), 1.62, "operating range", ha="center", fontsize=6, color="#444")
ax.text(2.2, 0.08, "stalled", fontsize=5.5, color=C_ARM)
ax.text(90, 1.25, "banked", fontsize=5.5, color=C_LEV)
ax.text(2.2, 0.75, "lever", fontsize=6.5, color=C_LEV); ax.text(200, 0.42, "spring-arm", fontsize=6.5, color=C_ARM)
ax.set_xscale("log"); ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
ax.set_xlabel("escape-wheel torque (µN·m)"); ax.set_ylabel("balance amplitude (rad)"); ax.set_ylim(-0.05, 1.8)
ax.set_title(f"Amplitude flat over {R['operating']['torque_dynamic_range']:.0f}× torque", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig2", "c")

# d: period error vs torque
fig, ax = new("fig2", "d")
for model, col in (("lever", C_LEV), ("arm", C_ARM)):
    ok = sw[(sw.model == model) & (sw.status == "ok")]
    ax.plot(ok.T_e*1e6, ok.rate_error*100, "-o", ms=2.5, color=col)
ax.axvspan(lo_T, hi_T, color="#e4e4e4", lw=0); ax.axhline(0, color=META_GREY, lw=0.6)
ax.text(3.0, 0.14, "lever", fontsize=6.5, color=C_LEV); ax.text(120, 0.075, "spring-arm", fontsize=6.5, color=C_ARM)
ax.set_xscale("log"); ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
ax.set_xlabel("escape-wheel torque (µN·m)"); ax.set_ylabel("period error (%)"); ax.set_ylim(-0.01, 0.19)
ax.set_title("Rate error < 0.1 % in range", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig2", "d")

# e: vs Q
fig, axs = plt.subplots(2, 1, figsize=(PX["fig2"]["e"][0]/DPI, PX["fig2"]["e"][1]/DPI), dpi=DPI, sharex=True)
axs[0].plot(qt.Q, qt.amplitude, "-o", ms=2.5, color=C_ARM); axs[0].set_ylabel("amplitude\n(rad)"); axs[0].set_ylim(0, 1.3)
axs[0].axhline(dp.THETA_MAX, color=META_GREY, lw=0.6, ls="--"); axs[0].text(42, dp.THETA_MAX+0.05, "banking", fontsize=5.5, color=META_GREY)
axs[1].plot(qt.Q, qt.rate_error*100, "-o", ms=2.5, color=C_ARM); axs[1].set_ylabel("period\nerror (%)"); axs[1].set_ylim(0, 0.5)
axs[1].set_xscale("log"); axs[1].set_xticks([40, 100, 400]); axs[1].set_xticklabels(["40", "100", "400"]); axs[1].minorticks_off(); axs[0].minorticks_off(); axs[1].set_xlabel("balance quality factor Q")
for a_ in axs: a_.yaxis.set_label_coords(-0.22, 0.5)
axs[0].set_title("Tolerant to unknown damping", loc="left")
fig.tight_layout(pad=0.3, h_pad=0.4); save(fig, "fig2", "e")

# ================================================================ FIG 3
# a: taut band schematic (side view)
fig, ax = new("fig3", "a"); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10.4); ax.set_aspect("equal"); fig.subplots_adjust(0, 0, 1, 1)
X = 6.3   # band axis
ax.add_patch(Rectangle((5.2, 9.4), 3.8, 0.5, fc=C_STEEL, ec="none")); ax.text(9.1, 9.65, "deck", fontsize=5.5, va="center")
ax.add_patch(Rectangle((3.6, 0.5), 5.4, 0.5, fc=C_STEEL, ec="none")); ax.text(9.1, 0.75, "base", fontsize=5.5, va="center")
ax.add_patch(Rectangle((X-0.45, 8.75), 0.9, 0.65, fc="#6b7680", ec="none")); ax.text(X+0.6, 9.05, "fixed clamp\n(slotted: L_h 18–26 mm)", fontsize=5.5, va="center")
ax.plot([X, X], [8.75, 5.95], color="#333", lw=1.5); ax.plot([X, X], [4.05, 1.95], color="#333", lw=1.5)
ax.text(X+0.25, 7.35, f"L_h = {dp.BAND_LH*1e3:.0f} mm\nT = T₀ + m g", fontsize=5.5, va="center"); ax.text(X+0.25, 3.0, f"L_h = {dp.BAND_LH*1e3:.0f} mm\nT = T₀", fontsize=5.5, va="center")
ax.add_patch(Rectangle((X-2.9, 4.05), 5.8, 1.9, fc="none", ec=C_PRINT, lw=1.8)); ax.add_patch(Rectangle((X-0.35, 4.05), 0.7, 1.9, fc=C_PRINT, ec="none"))
for xx in (X-2.4, X+2.4): ax.add_patch(Circle((xx, 5.0), 0.22, fc=C_BRASS, ec="none"))
ax.text(X-2.9, 6.25, "balance ring Ø62 mm, 7.4 g; M3 brass trim nuts ×4", fontsize=5.5, color=C_PRINT)
ax.add_patch(Rectangle((X-0.45, 1.3), 0.9, 0.65, fc="#6b7680", ec="none")); ax.text(X+0.6, 1.6, "tension clamp", fontsize=5.5, va="center")
ys = np.linspace(1.0, 1.3, 9); xs = X + 0.15*np.array([1,-1]*4+[1])
ax.plot(xs, ys, color="#555", lw=0.9); ax.text(X-0.6, 1.1, "coil spring T₀ = 1 N", fontsize=5.5, ha="right", va="center")
ax.text(0.1, 9.55, f"ribbon {dp.BAND_T*1e6:.1f} µm × 2 mm\n(0.0015\" feeler-gauge steel)", fontsize=6, va="top")
ax.text(0.1, 8.15, "k_half = GJ/L_h + T w²/(12 L_h)", fontsize=6, va="top")
ax.text(0.1, 7.45, f"k = {k*1e3:.2f} mN·m/rad\nf₀ = 1.20 Hz", fontsize=6, va="top")
ax.text(0.1, 3.6, f"no pivot, no lubricant;\ngravity only adds ±m g\nto one half-band\n→ ±{R['oscillator']['orientation_max_pct']:.2f} % in f", fontsize=6, va="top")
save(fig, "fig3", "a")
# b: frequency vs trim radius
rt = np.linspace(dp.R_TRIM_MIN, dp.R_TRIM_MAX, 200)
fs_ = np.array([ao.frequency(k, ao.balance_inertia(dp.BAL_R_OUT, dp.BAL_R_IN, dp.BAL_THICK, dp.RHO_PETG, dp.N_TRIM, dp.M_TRIM, r, dp.HUB_M, dp.HUB_R)[1]) for r in rt])
fig, ax = new("fig3", "b")
ax.axhspan(dp.F_TARGET*0.995, dp.F_TARGET*1.005, color=C_OK, lw=0)
ax.plot(rt*1e3, fs_, color=C_ARM, lw=1.8); ax.axhline(dp.F_TARGET, color=META_GREY, lw=0.7, ls="--")
ax.plot([r_reg*1e3], [dp.F_TARGET], "o", color=C_ALARM, ms=4); ax.text(r_reg*1e3+0.15, dp.F_TARGET+0.004, f"set point {r_reg*1e3:.1f} mm", fontsize=6, color=C_ALARM)
ax.text(dp.R_TRIM_MIN*1e3+0.1, dp.F_TARGET*1.0055, "±0.5 % calibration band", fontsize=6, color="#3a6b3a", va="bottom")
ax.text(dp.R_TRIM_MIN*1e3+0.1, dp.F_TARGET-0.0045, f"target {dp.F_TARGET:.4f} Hz (10 mL/h cassette)", fontsize=6, color=META_GREY, ha="left", va="top")
ax.set_xlabel("trim-nut radius (mm)"); ax.set_ylabel("balance frequency (Hz)")
ax.set_title(f"Trim weights regulate {R['oscillator']['regulation_range_pct'][0]:+.1f} % … {R['oscillator']['regulation_range_pct'][1]:+.1f} %", loc="left")
ax.margins(0.05); fig.tight_layout(pad=0.3); save(fig, "fig3", "b")

# c: orientation sensitivity lollipop on log axis
cs_sens, _, _ = acs.gravity_sensitivity(200e9, 5e-3, 0.05e-3, 10e-3, m_bal)
orient = R["oscillator"]["orientation_Hz"]; f_flat = orient["band horizontal (device flat)"]
vals = [("cross-spring flexure pivot\n(rejected design)", cs_sens/2*100, C_LEV),
        ("taut band, fixed clamp up", abs(orient["band vertical, fixed clamp up"]/f_flat-1)*100, C_ARM),
        ("taut band, fixed clamp down", abs(orient["band vertical, fixed clamp down"]/f_flat-1)*100, C_ARM)]
fig, ax = new("fig3", "c")
for i, (lab, v, col) in enumerate(vals):
    ax.plot([1e-3, v], [i, i], color=META_GREY, lw=0.8); ax.plot(v, i, "o", color=col, ms=5); ax.text(v*1.3, i, f"{v:.2g} %", va="center", fontsize=6)
ax.set_yticks(range(3)); ax.set_yticklabels([v[0] for v in vals]); ax.invert_yaxis()
ax.set_xscale("log"); ax.set_xlim(0.03, 40); ax.set_xticks([0.1, 1, 10]); ax.set_xticklabels(["0.1", "1", "10"])
ax.set_xlabel("frequency change vs. device flat (%)"); ax.set_title(f"{8.3/R['oscillator']['orientation_max_pct']:.0f}× less orientation error", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig3", "c")

# d: band stress vs amplitude
th = np.linspace(0, 1.4, 100)
fig, ax = new("fig3", "d")
ax.plot(th, ao.max_shear_stress(dp.BAND_T, th, dp.BAND_LH)/1e6, color=C_ARM, lw=1.8)
ax.axhline(ao.TAU_FATIGUE_STEEL/1e6, color=C_ALARM, lw=0.9); ax.text(0.02, ao.TAU_FATIGUE_STEEL/1e6+8, "fatigue limit, spring steel (fully reversed)", fontsize=6, color=C_ALARM)
ax.axvspan(qt.amplitude.min(), qt.amplitude.max(), color="#e4e4e4", lw=0); ax.text((qt.amplitude.min()+qt.amplitude.max())/2, 20, "operating\namplitude", ha="center", fontsize=6, color="#444")
ax.axvline(dp.THETA_MAX, color=META_GREY, lw=0.6, ls="--"); ax.text(dp.THETA_MAX+0.02, 60, "banking", fontsize=5.5, color=META_GREY)
ax.set_xlabel("balance amplitude (rad)"); ax.set_ylabel("peak shear stress (MPa)"); ax.set_ylim(0, 420)
ax.set_title("Infinite fatigue life", loc="left"); fig.tight_layout(pad=0.3); save(fig, "fig3", "d")

# e: temperature drift
Tc = np.linspace(10, 40, 50)
fig, ax = new("fig3", "e")
ax.plot(Tc, R["oscillator"]["tempco_per_K"]*(Tc-25)*100, color=C_ARM, lw=1.8)
ax.axhspan(-0.5, 0.5, color=C_OK, lw=0); ax.text(10.5, 0.42, "±0.5 %", fontsize=6, color="#3a6b3a", va="top")
ax.axhline(0, color=META_GREY, lw=0.6)
ax.text(10.5, -0.42, f"{R['oscillator']['tempco_per_K']*1e4:.2f} × 10⁻⁴ K⁻¹\n(steel G drift + ring expansion)", fontsize=6, va="bottom")
ax.set_xlabel("ambient temperature (°C)"); ax.set_ylabel("rate drift from 25 °C (%)"); ax.set_ylim(-0.6, 0.6)
ax.set_title("Drift ≤ 0.3 % over 15 K", loc="left"); fig.tight_layout(pad=0.3); save(fig, "fig3", "e")

# ================================================================ FIG 4
# a: flow heatmap
pv = ft.pivot(index="cassette", columns="syringe", values="flow_mL_h")
order_c = ["1 mL/h", "2 mL/h", "5 mL/h", "10 mL/h"]; order_s = ["BD 5 mL", "BD 10 mL", "BD 20 mL", "BD 30 mL", "BD 60 mL"]
pv = pv.loc[order_c, order_s]
fig, ax = new("fig4", "a")
im = ax.imshow(np.log10(pv.values), cmap="Blues", aspect="auto", vmin=np.log10(0.3), vmax=np.log10(25))
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]):
        v = pv.values[i, j]; ax.text(j, i, f"{v:.2g}", ha="center", va="center", fontsize=6.5, color="white" if v > 6 else "#222")
ax.set_xticks(range(5)); ax.set_xticklabels([s.replace("BD ", "") for s in order_s]); ax.set_yticks(range(4)); ax.set_yticklabels([c.replace(" mL/h", "") for c in order_c])
ax.set_xlabel("syringe (BD)"); ax.set_ylabel("cassette (mL/h on 20 mL)"); ax.set_title("Rates 0.4–19.5 mL/h (mL/h in cells)", loc="left")
for s in ax.spines.values(): s.set_visible(False)
fig.tight_layout(pad=0.3); save(fig, "fig4", "a")

# b: torque margin
fig, ax = new("fig4", "b")
cas = list(dp.CASSETTES); x = np.arange(len(cas))
Fmin, Fmax = opm.F_excess_N.min(), opm.F_excess_N.max()
for j, (sname, Pmm, col, lab) in enumerate((("BD 60 mL", 300.0, "#9ecae1", f"60 mL at 300 mmHg trip ({Fmin:.1f} N excess)"), ("BD 5 mL", 0.0, C_ARM, f"5 mL free-running ({Fmax:.1f} N excess)"))):
    d = opm[(opm.syringe == sname) & (opm.pressure_mmHg == Pmm)].set_index("cassette").loc[cas]
    ax.bar(x + (j-0.5)*0.38, d.margin, width=0.36, color=col, label=lab)
op1 = pd.read_csv("results/operating_map_rev1_30N.csv"); r1 = op1[(op1.syringe == "BD 60 mL") & (op1.pressure_mmHg == 300.0)].set_index("cassette").loc[cas]
ax.plot(x - 0.19, r1.margin, "x", color=C_ALARM, ms=5, label="rev 1, 30 N spring (stalled)")
ax.axhline(1, color=C_ALARM, lw=0.8); ax.text(-0.45, 0.78, "stall", color=C_ALARM, fontsize=6, ha="left", va="top")
ax.set_yscale("log"); ax.set_yticks([0.3, 1, 10, 100]); ax.set_yticklabels(["0.3×", "1×", "10×", "100×"]); ax.set_ylim(0.25, 3000)
ax.set_xticks(x); ax.set_xticklabels([c.replace(" mL/h", "") for c in cas]); ax.set_xlabel("cassette (mL/h)")
ax.set_ylabel("torque / re-cock minimum"); ax.set_title(f"Margin ≥ {R['operating']['worst_margin']:.1f}× in every cell", loc="left")
ax.legend(frameon=False, fontsize=5.2, loc="upper left", handlelength=1.2, borderaxespad=0.1); fig.tight_layout(pad=0.3); save(fig, "fig4", "b")

# c: resolution lollipop
res = R["resolution_uL_per_beat"]
fig, ax = new("fig4", "c")
for i, c in enumerate(cas):
    ax.plot([0, res[c]], [i, i], color=META_GREY, lw=0.8); ax.plot(res[c], i, "o", color=C_ARM, ms=5); ax.text(res[c]+0.05, i, f"{res[c]:.2f}", va="center", fontsize=6)
ax.set_yticks(range(len(cas))); ax.set_yticklabels([c.replace(" mL/h", "") for c in cas]); ax.invert_yaxis(); ax.set_xlim(0, 1.7)
ax.set_xlabel("µL per beat"); ax.set_ylabel("cassette (mL/h)"); ax.set_title("Sub-µL steps", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig4", "c")

# d: error budget
contrib = pd.Series(R["montecarlo"]["contrib_sd_pct"]).sort_values()
names = {"syringe": "syringe bore tolerance (±1 %)", "calibration": "trim regulation residual (±0.3 %)", "rack": "rack pitch after print calibration (±0.2 %)",
         "temperature": "ambient 15–35 °C", "orientation": "device orientation", "escapement": "escapement error, Q = 60–400"}
fig, ax = new("fig4", "d")
ax.barh([names[c] for c in contrib.index], contrib.values, color=[C_RESTR if c == "syringe" else C_ARM for c in contrib.index], height=0.6)
for y, v in enumerate(contrib.values): ax.text(v + 0.01, y, f"{v:.2f}", va="center", fontsize=6)
ax.set_xlabel("standard deviation of flow error (%)"); ax.set_title("Syringe, not mechanism, dominates", loc="left"); ax.set_xlim(0, 0.75)
fig.tight_layout(pad=0.3); save(fig, "fig4", "d")

# e: MC histogram (re-draw from the same generator as analysis_system)
rng = np.random.default_rng(7); N = 20000
d_cal = rng.uniform(-0.3, 0.3, N); Tamb = rng.uniform(15, 35, N); d_T = R["oscillator"]["tempco_per_K"]*(Tamb-25)*100
d_or = rng.uniform(-R["oscillator"]["orientation_max_pct"], R["oscillator"]["orientation_max_pct"], N)
Qs = rng.uniform(60, 400, N); esc = np.interp(np.log(Qs), np.log(qt.Q), qt.rate_error*100); d_esc = -(esc-esc.mean())
d_syr = rng.uniform(-1, 1, N); d_rack = rng.uniform(-0.2, 0.2, N); total = d_cal+d_T+d_or+d_esc+d_syr+d_rack
fig, ax = new("fig4", "e")
ax.axvspan(-5, 5, color=C_OK, lw=0, zorder=0); ax.hist(total, bins=80, color=C_ARM, alpha=0.9)
p2, p97 = np.percentile(total, [2.5, 97.5])
for p in (p2, p97): ax.axvline(p, color=META_GREY, ls="--", lw=0.7)
ax.set_ylim(0, ax.get_ylim()[1]*1.35)
ax.text(0.03, 0.95, f"95 % of units within\n{p2:+.2f} … {p97:+.2f} %\n(n = 20 000 draws)", transform=ax.transAxes, va="top", fontsize=6.5)
ax.text(4.9, ax.get_ylim()[1]*0.05, "±5 % requirement", ha="right", fontsize=6, color="#3a6b3a")
ax.set_xlim(-5.5, 5.5); ax.set_xlabel("steady-state flow-rate error (%)"); ax.set_ylabel("units"); ax.set_title("Predicted accuracy well inside ±5 %", loc="left")
fig.tight_layout(pad=0.3); save(fig, "fig4", "e")

# ================================================================ compose
paths = {}
for key, o in OUT.items():
    pp = {p["letter"]: f"panels/{key}_{p['letter']}.png" for p in o["panels"]}
    out, (Wp, Hp) = compose_figure(o, pp, f"manuscript/figures/{key}.png", letter_pt=11, letter_font=FONT_BOLD)
    paths[key] = out
print({k: v for k, v in paths.items()})
