
"""Fig 6: snap-action occlusion trip -- sensing travel, sear release, flag rod, stem capture; force vs travel."""
import json, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon, Circle
import design_params as dp
try:                      # figure helpers: kernel plugins (figure-style / figure-composer) or the standalone fallback
    apply_figure_style
except NameError:
    import os as _os, sys as _sys; _sys.path.append(_os.path.dirname(_os.path.abspath(__file__)) if "__file__" in globals() else _os.getcwd())
    from figstyle import apply_figure_style, panel_letter, META_GREY, end_of_line_labels, compose_figure, panel_px
R = json.load(open("results/system_results.json")); occ = R["occlusion"]["per_syringe"]
apply_figure_style(sizes=(8, 7, 6))
C_PRINT, C_STEEL, C_ALARM, C_ARM, META = "#1b6ca8", "#6b7680", "#c0392b", "#1b6ca8", "#7a7a7a"
fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1, 1])
ax = fig.add_subplot(gs[0, :]); ax.axis("off"); ax.set_xlim(0, 30.8); ax.set_ylim(0, 11.2); ax.set_aspect("equal")
def state(x0, title, sub, stem_dx, sear_up, flag_up, captured):
    ax.text(x0+0.3, 10.9, title, ha="left", fontsize=8, weight="bold", va="top"); ax.text(x0+0.3, 10.15, sub, ha="left", fontsize=5.8, va="top")
    # carriage body
    ax.add_patch(Rectangle((x0+0.3, 1.0), 3.2, 5.0, fc="#dbe4ea", ec=C_PRINT, lw=1.2))
    ax.text(x0+1.9, 1.3, "carriage", ha="center", fontsize=5.5, color=C_PRINT)
    # trip spring (zigzag) between carriage wall and stem collar; spacer
    sx0, sx1 = x0+3.5, x0+6.0 - stem_dx
    ax.add_patch(Rectangle((sx0, 3.3), 0.5, 1.4, fc="#e8a33d", ec="none")); ax.text(sx0+0.25, 2.95, "spacer", ha="center", fontsize=5, color="#b07a1a")
    xs = np.linspace(sx0+0.5, sx1, 13); ys = 4.0 + 0.45*np.array([0]+[1,-1]*5+[1,0])
    ax.plot(xs, ys, color=C_STEEL, lw=1.0)
    # stem + pusher plate
    ax.add_patch(Rectangle((sx1, 3.75), 2.2, 0.5, fc=C_STEEL, ec="none"))
    ax.add_patch(Rectangle((sx1+2.2, 2.4), 0.5, 3.2, fc=C_PRINT, ec="none")); ax.text(sx1+2.45, 2.0, "pusher", ha="center", fontsize=5.5, color=C_PRINT)
    # plunger reaction arrow
    ax.add_patch(FancyArrowPatch((sx1+3.6, 4.0), (sx1+2.75, 4.0), arrowstyle="-|>", mutation_scale=8, lw=1.2, color=C_ALARM))
    ax.text(sx1+3.2, 4.5, "F", fontsize=6, color=C_ALARM, ha="center")
    # sear (schematic): pawl above the stem; capture slot on the stem top (CAD: Fig. 6c,d)
    nx = sx1 + 0.9
    ax.add_patch(Polygon([[nx-0.25, 4.25], [nx+0.25, 4.25], [nx, 3.95]], fc="white", ec=C_STEEL, lw=0.6))  # notch
    sy = 5.0 + (0.6 if sear_up else 0.0)
    ax.add_patch(Polygon([[x0+6.6, 6.4], [x0+6.9, sy], [x0+7.3, sy], [x0+7.4, 6.4]], fc="#e8a33d", ec="#b07a1a", lw=0.8))  # sear
    ax.add_patch(Circle((x0+7.0, 6.4), 0.18, fc="white", ec="k", lw=0.6))
    if captured:
        ax.add_patch(Polygon([[x0+6.9, 4.25], [x0+7.3, 4.25], [x0+7.1, 3.95]], fc="#e8a33d", ec="#b07a1a", lw=0.6))
    # flag rod: vertical rod held by sear foot when armed; released -> up
    fy = 6.9 + (1.5 if flag_up else 0.0)
    ax.add_patch(Rectangle((x0+7.8, 5.2 + (1.5 if flag_up else 0)), 0.3, 1.8, fc=C_STEEL, ec="none"))
    ax.add_patch(Rectangle((x0+7.5, fy), 0.9, 0.6, fc=C_ALARM, ec="none"))
    # flag spring (compressed when armed)
    fs_y = np.linspace(5.2 + (1.5 if flag_up else 0), 4.3, 7); fs_x = x0+7.95 + 0.18*np.array([1,-1]*3+[1])
    ax.plot(fs_x, fs_y, color=C_STEEL, lw=0.8)
    # rail that blocks the balance
    ax.add_patch(Rectangle((x0+8.3, 5.5 + (1.5 if flag_up else 0)), 1.2, 0.25, fc=C_ALARM if flag_up else C_STEEL, ec="none"))
    ax.text(x0+9.0, 5.1 + (1.5 if flag_up else 0), "rail", fontsize=5, ha="center", color=C_ALARM if flag_up else META)
    # balance ring edge
    ax.add_patch(Circle((x0+9.9, 8.0), 0.8, fc="none", ec=C_ARM, lw=1.3)); ax.text(x0+9.9, 6.6, "balance rim", fontsize=5, ha="center", color=C_ARM)
state(0.0, "1  armed", "preload $F_0$ = spacer × 3 N/mm;\nsear holds flag rod down", 0.0, False, False, False)
_ds, _dF, _fl = dp.TRIP_SENSE_TRAVEL, dp.TRIP_SPRING_RATE*dp.TRIP_SENSE_TRAVEL, dp.FLAG_TRAVEL
state(10.3, "2  sensing", f"F > $F_0$: stem retracts {_ds:.2f} mm\n(+{_dF:.2f} N); trip pin turns the sear", _ds*1.5, False, False, False)
state(20.6, "3  tripped", f"sear frees flag rod (up {_fl:.0f} mm): rocker drops capture\npin into stem slot; rail lifts fork onto balance spoke", _ds*1.5, True, True, True)
ax.text(0.2, 0.3, "stem captured ⇒ stored spring compression cannot advance the plunger when the blockage clears; reset by the ramp at the loading (fully retracted) position", fontsize=5.5, color="#333333")
# panel b: force vs stem travel
ax2 = fig.add_subplot(gs[1, 0])
k = R["occlusion"]["spring_rate_N_per_mm"]; ds = R["occlusion"]["sense_travel_mm"]
x = np.linspace(0, 3.2, 200)
for sname, col in (("BD 20 mL", C_ARM), ("BD 60 mL", "#2c7fb8"), ("BD 5 mL", "#9ecae1")):
    F0 = occ[sname]["spacer_preload_N"]; Ft = occ[sname]["trip_force_N"]
    ax2.plot(x, F0 + k*x, color=col, lw=1.4)
    ax2.plot([ds], [Ft], "o", color=col, ms=5)
    ax2.text(3.25, F0 + k*3.2, sname.replace("BD ", ""), fontsize=6, color=col, va="center")
ax2.axvline(ds, color=C_ALARM, lw=0.9); ax2.text(ds+0.08, 33, f"sear releases\n({ds:.2f} mm)", fontsize=6, color=C_ALARM, va="top")
ax2.axvspan(ds, 3.0, color="#eeeeee", lw=0, zorder=0); ax2.text(1.7, 2.5, "rev 2: 3 mm direct\ntravel (+9 N; +235 mmHg\non 20 mL)", fontsize=5.5, color=META, ha="center")
ax2.set_xlabel("stem travel (mm)"); ax2.set_ylabel("trip-spring force (N)"); ax2.set_xlim(0, 3.9); ax2.set_ylim(0, 37)
ax2.set_title(f"Trip force is defined at {ds:.2f} mm, not 3 mm")
# panel c: per-syringe trip table
ax3 = fig.add_subplot(gs[1, 1]); ax3.axis("off")
rows = [["syringe", "preload\n(N)", "trip\n(N)", "fires at\n(mmHg)", "stored at sear\n(mL)"]]
for sname, v in occ.items():
    lo, nom, hi = v["pressure_band_mmHg"]
    rows.append([sname.replace("BD ", ""), f"{v['spacer_preload_N']:.1f}", f"{v['trip_force_N']:.1f}", f"{lo:.0f}–{hi:.0f}", f"{v['stored_displacement_at_sear_mL']:.2f}"])
tb = ax3.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center", colLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(6); tb.scale(1.0, 1.45)
for c_ in range(5): tb.auto_set_column_width(c_)
for (r_, c_), cell in tb.get_celld().items():
    cell.set_edgecolor("#bbbbbb"); cell.set_linewidth(0.5)
    if r_ == 0: cell.set_text_props(weight="bold"); cell.set_height(cell.get_height()*1.6)
ax3.set_title("Trip settings per syringe (300 mmHg nominal)")
ax3.text(0.5, -0.02, "fires-at band: ±1.5 N plunger-friction uncertainty.\nIdeally the captured stem releases none of the stored displacement\n(geometric restraint verified by solid Booleans, Fig. 7d); the residual\ndisplacement and total release bolus are measured in protocol step 5", transform=ax3.transAxes, ha="center", va="top", fontsize=6, color="#444444")
panel_letter(ax, "a"); panel_letter(ax2, "b"); panel_letter(ax3, "c")
fig.savefig("figures/fig6_occlusion_trip.png", dpi=300)
