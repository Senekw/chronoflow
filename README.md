# ChronoFlow

**An electricity-free, volumetric, clockwork-governed syringe infusion driver for neonatal care in low-resource settings.**

Design-and-simulation study: parametric 3D CAD (95 parts, interference-checked in 12 poses), physics-based validation, LaTeX manuscript, ISEF abstract and a bench validation protocol. *No physical prototype has been built yet.*

![overview](figures/manuscript/fig1.png)

## What it is
A wind-up syringe driver that meters **volume**, not pressure. Two constant-force springs (44.5 N) push the plunger; a rack-and-pinion drives a printed 312.5:1 step-up train into a clockwork escapement whose balance wheel sets the rate. Plug-in gear cassettes give 1 / 2 / 5 / 10 mL/h on a BD 20 mL syringe (0.4–19.5 mL/h across 5–60 mL syringes).

Two design results carry the project:
1. **Constant-force spring-arm escapement** — train torque varies 24× across syringe sizes, rates and loads (19–476 µN·m); a conventional lever escapement over-banks, the spring-arm design holds amplitude constant with a period error of 0.06 %.
2. **Taut-band balance** — a 38.1 µm (0.0015") × 2 mm feeler-gauge steel ribbon under 1 N pretension replaces pivots: 1.20 Hz, ±0.19 % orientation sensitivity (a cross-spring flexure pivot would be ~8 %), 166 MPa vs 350 MPa fatigue limit; a slotted clamp (L_h 18–26.5 mm) absorbs the ±2.5 µm stock tolerance.

Predicted accuracy: 95 % of units within −1.14 … +1.15 % (requirement ±5 %); flow-restriction pumps: +45 % at 37 °C, −33 % at 100 mmHg back-pressure. Escapement torque margin ≥ 2.0× over all 5 syringes × 4 cassettes × 0–300 mmHg. BOM ≈ $40.

**Revision 5** (external review round 6): arm-to-detent release linkage (release tails + release star + radial detent lever) and the trip reset ramp are now modelled and Boolean-checked; impulse energy predicted forward from a specified wire torsion spring; printed arm flexures replaced (they are 150× too stiff, not too soft); re-cock budget restated as half period minus impulse (0.36 s, margin 2.5×); "zero release bolus by construction" withdrawn in favour of a measured residual bolus; J_eff, detent deflection, protocol and bore-tolerance statements corrected.

**Revision 4** (external review round 5): the snap-action occlusion trip is now fully modelled in CAD (keyed spacer, trip pin, sear on a vertical pivot, 2 N flag rod, rocker-driven capture pin, trip rail on parallel-motion cranks, drop rod and balance stop fork) in an *armed* and a *tripped* pose, and the stem capture is verified on the solids; the band clamps are slotted (hanger plate + bracket, M3 screws); a mechanism-level escapement analysis (`analysis_arm_mechanism.py`, Fig. 3) derives impulse and unlatch energies from geometry, integrates the wheel through re-cocking and sizes a compliant wheel detent (now in CAD); the interference check uses an explicit pair whitelist over 12 poses with a carriage-posed gear train; all 91 STLs are watertight; the figure scripts run standalone (`figstyle.py`). See `docs/ChronoFlow_Report.md` for the itemised log.

**Revision 2** (external review): drive spring 30 → 44.5 N (the 60 mL syringe stalled before its occlusion trip at 30 N), snap-action occlusion trip (0.3 mm sensing travel to a sear, stem captured, 12 mm sprung flag rod) with spacer preload keyed to the syringe cradle insert and a disclosed friction-dependent pressure band, both band clamps slotted for symmetric coarse regulation, ribbon modelled at the true 38.1 µm stock thickness with a 22 mm half-band, mesh count corrected (rack mesh at 0.95 plus six printed gear meshes at 0.92, not seven at 0.92). See `docs/ChronoFlow_Report.md` §Revision 2.

## Repository layout
| Path | Contents |
|---|---|
| `analysis/design_params.py` | single source of truth for every dimension, ratio and spring |
| `analysis/analysis_oscillator.py` | taut-band stiffness, frequency, orientation, stress, temperature |
| `analysis/analysis_crossspring.py` | beam-column FE of the rejected cross-spring pivot |
| `analysis/analysis_escapement.py` | RK4 balance–escapement simulator (lever vs spring-arm) |
| `analysis/analysis_system.py` | flow tables, sweeps, operating map, Monte-Carlo, BOM → `results/`, analysis figures |
| `analysis/analysis_arm_mechanism.py` | mechanism-level escapement: impulse/unlatch energies from geometry, wheel re-cock ODE, detent impact, rate sensitivity |
| `analysis/make_panels.py`, `make_fig6.py` | manuscript figure panels + composition |
| `analysis/figstyle.py` | standalone copy of the plotting helpers (`apply_figure_style`, `panel_letter`, `compose_figure`, …) used when the kernel plugins are absent |
| `cad/cad_model.py` | CadQuery parametric model → STEP / glTF / STL |
| `cad/check_interference.py` | pairwise OCC Boolean interference check: explicit allowed-pair list, 12 poses, trip-capture functional check → `results/interference_report.json`, `results/trip_capture_check.json` |
| `cad/render_cad.py`, `cad/make_renders.py` | OpenGL-free painter's renderer and the script that renders the manuscript CAD panels (armed + tripped manifests) |
| `cad/chronoflow_assembly.step`, `.glb`, `cad/stl/`, `cad/stl_tripped/` | exported geometry (all STLs watertight; tripped pose exported for the figure) |
| `viewer/chronoflow_viewer.html` | self-contained three.js viewer (orbit, click-to-identify, hide deck, explode) |
| `manuscript/chronoflow.tex` | full manuscript (figures in `manuscript/figures/`) |
| `docs/` | long-form design report and ISEF-format abstract |
| `results/` | CSV/JSON outputs behind every number |

## Reproduce
```bash
conda create -n chronoflow python=3.11 numpy scipy pandas matplotlib trimesh -c conda-forge
conda activate chronoflow && pip install cadquery
cd analysis && python analysis_system.py          # results/*.csv|json + figures/analysis/*.png
python analysis_arm_mechanism.py                  # results/arm_mechanism.json + fig7
python make_fig6.py && python make_panels.py      # manuscript figure panels + composed figures
cd ../cad && python cad_model.py                  # STEP + glTF + STLs (armed and tripped) + manifests (≈15 s)
python check_interference.py                      # 12-pose interference check + trip-capture check (≈2.5 min)
python make_renders.py                            # CAD panels for Fig. 1a and Fig. 7a–e
```
The figure scripts import the plotting helpers from `analysis/figstyle.py` when the figure-style / figure-composer kernel plugins are not loaded, so they run from a plain shell. Run the scripts from their own directory (they resolve `results/`, `figures/` and `panels/` relative to the repository root as laid out here).

## Key numbers (from `results/*.json`, regenerated with this revision)
- Oscillator: k = 0.295 mN·m/rad, I = 5.18e-06 kg·m², f = 1.2017 Hz, trims at 26.4 mm, regulation -2.4 … +3.1 %, orientation ±0.19 %, band stress 166 MPa
- Escapement: T_min = 9.5 µN·m, worst margin 2.04× (BD 60 mL, 1 mL/h, 300 mmHg), torque 19.5–476 µN·m, power 17–416 µW; rate error 0.057 % at Q = 150
- Mechanism level: impulse energy predicted forward from the spring (k = 1e-4 N·m/rad, 26.3° preload, µ = 0.15): E_imp = 1.20 µJ, droop 6 %; re-cock + lock ≤ 0.14 s vs the 0.36 s left in a beat interval (half period − impulse; worst 1 mL/h, margin 2.53×, 3.2–4.0× other cassettes); unlatch 0.10–0.35 µJ (µ 0.1–0.4); surplus ≤ 98 µJ/beat → detent impact 0.62 N / 19 MPa / 0.31 mm with the flexure detent (14 N / 92 MPa if rigid); detent-release pin withdrawal 0.94 / 0.83 mm vs 0.6 mm required (margin 1.6 / 1.4); period error 0.050–0.108 % over E_imp ±10 % and the unlatch range
- Occlusion trip: sensing travel 0.35 mm (sear geometry), preloads 5.0–25.8 N, flag travel 6 mm; capture check: pin engagement 1.45 mm, stem advanced 0.1 mm → 0.085 mm³ pin/stem interference (captured); balance free at ±0.04 rad, blocked at ±0.08 rad inside the fork; reset ramp at the loading position: tripped head/ramp overlap 117 mm³, armed head clear
- Structure: rack-pinion tooth stress 22.7 MPa (SF 2.2)
- Energy: 4.23 J per wind, 0.37 J used in 20 h
- CAD: 95 parts, 0 interference faults in 12 poses (28 allowed pairs with reasons), 95/95 STLs watertight

## Status and limitations
All performance figures are model predictions. The balance quality factor in air and printed-gear efficiencies are the least certain inputs (both swept). The spring-arm escapement in the CAD is a representative embodiment that will need bench iteration; the detent-release margin (1.4–1.6×) and the impulse contact parameters (preload, µ) are the first things to measure. See `docs/ChronoFlow_Report.md` §5–6 for the validation protocol and limitations.

## License
MIT (code) · CC-BY-4.0 (CAD, figures, text). See `LICENSE`.
