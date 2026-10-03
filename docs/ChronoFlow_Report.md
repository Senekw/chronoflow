# ChronoFlow: an electricity-free, volumetric, clockwork-governed syringe infusion driver for neonatal care in low-resource settings

**Category:** Engineering — Biomedical Engineering / Engineering Mechanics
**Status:** Complete design, parametric 3D CAD and physics-based validation. No physical prototype has been built yet; a build-and-test protocol is provided.

---

## Abstract

Hospitals without reliable electricity deliver intravenous drugs and fluids to newborns by gravity drip or by hand, with dosing errors that are often fatal at neonatal volumes. Existing non-electric infusion devices meter flow by pushing fluid through a narrow restriction, so the delivered rate depends on fluid viscosity (temperature, drug) and on downstream pressure; electric low-cost pumps still need batteries and cost hundreds of dollars. ChronoFlow is a wind-up syringe driver that meters *volume*, not pressure: a constant-force spring pushes the plunger, and a rack-and-pinion drives a 3D-printed step-up gear train into a clockwork escapement whose balance wheel sets the rate. Two design problems were identified and solved by simulation. (1) The mechanical power available to the escapement varies ~24× across syringe sizes, clinical rates and loads, so a conventional lever escapement over-banks or stalls; adapting the 19th-century gravity-escapement principle with springs (a *constant-force spring-arm escapement*) delivers a fixed impulse energy per beat and holds the balance amplitude constant over the whole range (period error < 0.1 %). (2) Pivot friction and orientation sensitivity were eliminated with a *taut-band torsion suspension* made from a 0.0015" (38.1 µm) feeler-gauge steel ribbon (frequency shift ±0.19 % between upright and inverted). Interchangeable colour-coded gear cassettes give 1, 2, 5 and 10 mL/h with a 20 mL syringe (0.4–19.5 mL/h across 5–60 mL syringes) at 0.12–1.16 µL per beat resolution. A 20 000-draw Monte-Carlo error budget predicts 95 % of units within ±1.14 % of set rate (requirement ±5 %), compared with +45 % at 37 °C and −33 % at 100 mmHg back-pressure for flow-restriction pumps. A trip spring whose preload is keyed to the syringe cradle insert stops the oscillator and raises a flag at 300 mmHg occlusion (280–320 to 202–398 mmHg across syringe sizes once plunger-friction uncertainty is included), and a syringe-latch interlock prevents free-flow boluses. The complete 66-part assembly was modelled parametrically (CadQuery), checked for interferences, and exported as STEP/STL; estimated bill of materials is $40.

---

## 1. Problem and prior art

Neonates need fluids at 60–100 mL/kg/day and drugs such as antibiotics, dextrose and inotropes at 0.5–10 mL/h; a 10 % rate error or an uncontrolled bolus of a few millilitres is clinically significant. The device landscape divides into:

| Approach | Example | Limitation |
|---|---|---|
| Electronic syringe pump | commercial | $1 000–3 000, mains/battery, service infrastructure |
| Low-cost electro-mechanical | AutoSyP (Rice 360), constant-force spring + motor | ≈$500, rechargeable battery, 3–60 mL/h |
| Flow-restriction, non-electric | Springfusor, elastomeric pumps, air-over-oil design (ASME J. Med. Dev. 2017) | rate ∝ ΔP/μ: depends on viscosity, temperature, back-pressure; proprietary consumables |
| Gravity + pendulum metronome | ASME 2011 student design | orientation-dependent, needs a still, level surface |
| Kitchen-timer clockwork | microfluidics, 2024 | one fixed rate, ≈1 h runtime |

The gap: a non-electric device that meters **volume**, works in **any orientation**, offers **selectable clinical rates**, runs for a full syringe without re-winding, and has **mechanical occlusion detection and free-flow protection**.

## 2. Design

![Figure 1 — labelled overview]({{artifact:art_c6864693-93a7-43e4-9a3c-2a3bcd6c0ec5}})

*Figure 1. ChronoFlow assembly (BD 20 mL syringe at 85 % fill, 10 mL/h cassette). Upper deck: syringe cradle with flange latch and barrel clamp, carriage on two 6 mm rails with an integrated module-1 rack, pusher head with occlusion trip spring and flag, constant-force spring drum. The rack pinion carries a one-way click so that pulling the carriage back to load a syringe is also the winding action. Below the deck: 5 × 5 × 5 × 2.5 step-up train, rate cassette, escape wheel, spring-arm escapement and the taut-band balance (visible through a side window so staff can see that it is running).*

### 2.1 Architecture and energy flow

Two constant-force springs (2 × 5 lbf, 44.5 N total) pull the carriage toward the syringe nozzle. The plunger load (syringe friction 1.5–4.5 N plus pressure × bore area, up to 26.9 N for a 60 mL syringe at its 300 mmHg trip) absorbs part of this force; the *excess* back-drives the rack pinion through the step-up train to the escape wheel. The escapement releases one half-tooth per beat of the balance, so plunger displacement per beat is fixed by geometry:

Δx per beat = π·m·z₀ / (2·Z_e · N_fixed · r_cassette) = π·1·12 / (2·15·312.5·1) = 4.0 µm (10 mL/h cassette), i.e. 1.16 µL per beat on a 20 mL syringe (0.12 µL at 1 mL/h).

Because the metered quantity is displacement, fluid viscosity, temperature and downstream pressure do not enter the rate (Fig. 3) until the occlusion trip point. The spring never "runs down": force is constant over the whole 95 mm travel, so runtime is set only by syringe volume (2–20 h per 20 mL fill).

### 2.2 Oscillator: taut-band balance

A 62 mm printed PETG ring (7.4 g, I = 5.2 × 10⁻⁶ kg m²) is clamped at the midpoint of a 38.1 µm × 2 mm spring-steel ribbon (a 0.0015" feeler-gauge blade; revision 1 modelled 40 µm, 6 % off in frequency) held at both ends and pre-tensioned to 1 N by a steel coil spring. Each 22 mm half-band is a torsion spring, k_half = GJ/L + T w²/(12 L) (Timoshenko), giving k = 2.95 × 10⁻⁴ N m rad⁻¹ and f₀ = 1.20 Hz. Because k ∝ t³, the ±2.5 µm tolerance of feeler stock alone moves f₀ by ±9 %, so regulation is two-stage: both end clamps sit on slotted posts and are moved symmetrically (hub stays centred; k = k_half(L_up) + k_half(L_low) with each half's actual length), so that L_h = 18.4–26.1 mm re-centres f₀ for any ribbon between 35.6 and 40.6 µm once the ribbon has been measured — adjusting one clamp alone halves the sensitivity (1.1 %/mm vs 2.2 %/mm) and the same slot would cover only 36.6–39.4 µm — and four brass M3 nuts in radial slots (23–29 mm) then regulate by −2.4 % … +3.1 %, exactly as a clock is regulated. Design checks:

| Property | Value | Note |
|---|---|---|
| Frequency | 1.2017 Hz (target 1.2017) | L_h = 22 mm, trims at 26.4 mm |
| Orientation sensitivity | ±0.19 % | gravity changes tension of one half-band only |
| Peak shear stress at 1.2 rad | 166 MPa | fatigue limit ≈ 350 MPa → infinite life |
| Temperature coefficient | −1.8 × 10⁻⁴ K⁻¹ (−0.45 % over 25 K) | steel G drift + ring expansion |
| Band sag (device on its side) | 0.40 mm | within clearances |

An earlier candidate, a cross-spring flexure pivot, was rejected by analysis: a finite-element beam-column model showed that gravity on a 10 g balance changes a 0.05 mm leaf pivot's rotational stiffness by ≈20 % (the relevant scale is EI/L² ≈ 0.1 N, not the Euler load), which would make the rate orientation-dependent.

### 2.3 Escapement: constant-force spring-arm

The escapement must run from 17 µW (60 mL syringe, 1 mL/h cassette, at its 300 mmHg occlusion trip) to 416 µW (5 mL syringe, 10 mL/h, empty line) — a 24× torque range at the escape wheel (19–476 µN m). A balance driven by a conventional lever escapement receives an impulse proportional to that torque, so its amplitude scales as √torque: the simulation (Fig. 2a) shows it banking at > 1.2 rad over most of the operating range, which both damages rate and risks over-banking. The ChronoFlow escapement instead uses two spring-loaded impulse arms. Each escape-wheel tooth *cocks* an arm against a latch; when the balance passes, its impulse pin trips the latch and the arm delivers a fixed energy E_imp = 1.2 µJ over a 15° arc. Train torque only re-cocks the arm, so any excess is dissipated in the locking impact and never reaches the balance — the spring analogue of Denison's gravity escapement in tower clocks, whose purpose was exactly to isolate the pendulum from variable train torque.

![Figure 2 — escapement comparison]({{artifact:art_09f770dc-3d9a-4622-802c-b8ca2ef980ed}})

*Figure 2. Balance–escapement dynamics (RK4, Q = 150). (a) Amplitude vs escape-wheel torque: lever escapement amplitude grows as √T and exceeds the banking limit inside the operating range; the spring-arm escapement holds 0.59 rad from the minimum re-cock torque upward. (b) Period error relative to the free oscillator stays below 0.06 % for the spring-arm design (Airy escapement error from the asymmetric impulse window).*

Minimum re-cock torque is T_min = E_imp/(η_lift·π/Z_e) = 9.5 µN m; the carriage speed is set by the cassette alone, so torque depends only on cassette and excess force, and the operating map was evaluated over all 5 syringes × 4 cassettes × back-pressures 0–300 mmHg (80 cells). The worst cell (60 mL syringe, 1 mL/h cassette, at its 26.9 N trip) supplies 19.5 µN m, a 2.0× margin (Fig. 5). This map is what raised the drive spring from the 30 N of revision 1 to 44.5 N: at 30 N the 60 mL trip left 3.1 N excess, 3.4 µN m at the slowest cassette, and the escapement would have stalled before the occlusion flag was raised (two cells stalled, four below 1.7×). With unknown air damping, the amplitude settles between 0.31 rad (Q = 40) and 0.99 rad (Q = 800) while the period error stays within 0.01–0.45 %, which is removed by trim calibration.

### 2.4 Rate selection: gear cassettes

A plug-in cassette carries two gear pairs with constant tooth sums (z₁+z₂ = z₃+z₄ = 66, 33 mm centres): (33,33,33,33) → 10 mL/h, (44,22,33,33) → 5, (33,33,55,11) → 2, (44,22,55,11) → 1 mL/h on a 20 mL BD syringe. The cassette is a visible, colour-coded physical object, so the set rate can be confirmed at a glance and cannot drift. Flow scales with syringe bore for other sizes:

| Cassette | BD 5 mL | BD 10 mL | BD 20 mL | BD 30 mL | BD 60 mL |
|---|---|---|---|---|---|
| "10" | 3.97 | 5.75 | 10.0 | 12.9 | 19.5 |
| "5" | 1.99 | 2.87 | 5.0 | 6.4 | 9.7 |
| "2" | 0.79 | 1.15 | 2.0 | 2.6 | 3.9 |
| "1" | 0.40 | 0.57 | 1.0 | 1.3 | 1.95 |

(mL/h; the 11-tooth pinion should be cut with +0.3 profile shift.)

### 2.5 Safety functions (all mechanical)

* **Occlusion trip (snap-action, Fig. 6; CAD in Fig. 7c–d).** The pusher head rides on an Ø8 steel stem backed by a 3 N/mm compression spring whose preload is set by a spacer keyed to each syringe cradle insert. When plunger reaction exceeds the preload the stem retracts; a trip pin on its end collar stands 0.15 mm behind the tail of a sear (a printed plate on a vertical pivot whose hook lip overlaps the flag-rod collar by 0.9 mm; lever ratio 4.5), so the *sensing travel* is 0.15 + 0.9/4.5 = 0.35 mm — a property of the modelled parts, not an assumption. The freed 2 N flag rod rises 6 mm: it exposes the red flag, releases a torsion-sprung rocker that drives a Ø2 capture pin 2.2 mm into a slot in the stem top, and lifts a frame-mounted rail (on two parallel-motion cranks, so it translates 6 mm wherever the carriage is) whose drop rod raises a two-pin fork around a balance spoke, confining the balance to ±0.05 rad (< the 6° unlock angle): ticking stops and pressure cannot rise further. Defining the trip at 0.35 mm rather than 3 mm of direct travel matters: 3 mm would add 9 N (+235 mmHg on 20 mL, and enough to stall the 60 mL case before the flag appeared); 0.35 mm adds 1.05 N, pre-compensated in the spacer so the sear releases at 300 mmHg × A + friction — trip forces 6.1/8.6/14.5/18.3/26.9 N from spacer preloads 5.0/7.6/13.4/17.2/25.8 N (1.7–8.6 mm). The slot's forward wall is cut 0.35 mm short of the pin so that after the sensing travel it is flush with it: the 0.04–0.20 mL stored in the spring cannot be released as a bolus. This is verified on the solids (`results/trip_capture_check.json`): advancing the stem 0.1/0.2 mm in the tripped pose gives 0.085/0.235 mm³ pin–stem interference, the armed and tripped poses are interference-free, and the balance is free at ±0.04 rad and blocked at ±0.08 rad inside the fork. For ±1.5 N plunger-friction uncertainty the 300 mmHg setting fires at 280–320 mmHg (60 mL) to 202–398 mmHg (5 mL); the protocol measures glide force per lot.
* **Free-flow protection.** The LOAD/RUN lever opens the flange latch *and* lifts the rack-pinion click at the same time. The carriage can be moved by hand only when the syringe is unlatched, so a fast-forwarded pusher moves the whole syringe along the cradle rather than expelling fluid.
* **Winding = loading.** Pulling the carriage back against the 44.5 N spring both loads a fresh syringe and stores 4.2 J — about eleven times the 0.37 J the escapement consumes in a 20 h infusion.
* **Audible/visual running indication.** The balance beats at 1.2 Hz behind a window; silence means stopped.

### 2.6 Structural check

The highest-loaded teeth are the rack and 12-tooth rack pinion: Lewis bending stress at 44.5 N on an 8 mm face is 22.7 MPa against ≈50 MPa PETG yield (safety factor 2.2; a 10 mm stage-1 face would restore 2.75). Torque falls 5× per stage (267 → 49 → 9.0 → 1.7 mN m), so the remaining printed gears are lightly loaded; shafts are 3 mm steel pins in printed bores.

## 3. Predicted performance

![Figure 3 — comparison]({{artifact:art_028178a8-f664-4eff-9b94-4ec5011cb8a3}})

*Figure 3. (a) Flow-restriction pumps follow Q ∝ ΔP/μ: water viscosity (Vogel equation) falls 31 % from 20 °C to 37 °C, so a pump calibrated at 20 °C over-delivers by 45 % at body temperature and under-delivers 23 % at 10 °C; ChronoFlow drifts −0.36 % … +0.18 %. (b) A 100 mmHg venous/catheter back-pressure cuts a 300 mmHg restriction pump by 33 %; ChronoFlow is unaffected until it trips at 300 mmHg.*

![Figure 4 — error budget]({{artifact:art_cad7eac5-9c28-41b0-b52a-c9bcd2044302}})

*Figure 4. Monte-Carlo steady-state flow error (n = 20 000). Inputs: syringe bore ±1 % (uniform), trim regulation residual ±0.3 %, rack pitch after print-scale calibration ±0.2 %, ambient 15–35 °C, orientation ±0.19 %, escapement error for Q = 60–400. Result: 95 % interval −1.14 … +1.15 %, 99th percentile |error| 1.36 %. The dominant term, syringe bore tolerance, is common to every syringe pump including electronic ones.*

![Figure 5 — operating map]({{artifact:art_f7555092-cefb-404b-8274-e5a3aa1d1b45}})

*Figure 6. Snap-action occlusion trip: armed / sensing / tripped states, trip-spring force vs stem travel (sear at 0.35 mm; previous 3 mm direct-travel design shaded), and per-syringe spacer preload, trip force, firing band and displacement stored at the sear.*

*Figure 5. (a) Escape-wheel torque relative to the arm re-cock minimum for every syringe × cassette combination at the 300 mmHg trip. (b) The two extremes of the map — 60 mL at trip (17.6 N excess) and 5 mL free-running (43.0 N excess) — with the first-revision 30 N spring shown as crosses (stalled at 2 and 1 mL/h). Train efficiency 0.95 × 0.92⁶ = 0.58: seven meshes in total — the rack mesh at 0.95 and six printed gear meshes at 0.92 (four step-up, two cassette); the escape wheel shares the cassette output shaft, so there is no seventh gear mesh.*

Resolution: 1.16 / 0.58 / 0.23 / 0.12 µL per beat (10/5/2/1 mL/h), i.e. 2.4 micro-boluses per second — smoother than a gravity drip (≈50 µL drops) and comparable to stepper-driven pumps.

## 4. CAD model and manufacturability

![Escapement close-up]({{artifact:art_e9ee872b-977c-49e5-878b-42357f2c15b1}})

![Exploded view]({{artifact:art_f029c561-96b4-44c1-b30f-ee6b681db7e4}})

The model is fully parametric: `design_params.py` holds every dimension and ratio; `cad_model.py` (CadQuery) builds the 91 parts in a common frame, derives shaft positions from the gear centre distances, phases every gear pair so a wheel tooth faces a pinion gap, poses the whole train with the carriage position, builds an *armed* and a *tripped* state of the occlusion mechanism, and exports STEP, glTF and per-part STL (all 91 watertight). `check_interference.py` runs a pairwise solid-interference check (bounding-box prefilter + OCC Boolean common, 4 095 pairs) in 12 poses — carriage full/half/empty × armed/tripped × balance at rest and ±0.6 rad — against an explicit list of 20 allowed part pairs (press-fit pins, threaded screws, the drawn 0.1 mm penetration of the locked tooth, the ≤ 0.7 mm³ tip overlap of the unshifted 12T pinions), each with a reason and a volume cap; the last run reports zero faults (`results/interference_report.json`). The 12T pinion tip overlap is the known undercut of a 12-tooth 20° pinion (remedy: +0.3 profile shift). All printed parts are FDM-friendly (module 1 teeth, ≥ 2 mm walls, 0.8 mm flexures) and print in ≈420 g of PETG. Bill of materials ≈ $38 (springs $13, printed parts $9, rails $3.60, cassettes $4.80, hardware/pins/feeler stock $7.80).

Interactive model: `chronoflow_viewer.html` (orbit, part identification, hide deck, explode).

## 5. Build-and-test protocol (next phase)

1. **Oscillator bench.** Print the balance; clamp the ribbon; measure f₀ with a phone slow-motion camera (240 fps) or a light gate. Target 1.20 Hz; set trims. Measure Q from the free-decay envelope (expected 60–400). Repeat with the band horizontal, fixed-clamp-up and fixed-clamp-down; acceptance ±0.3 %.
2. **Escapement bench.** Drive the escape wheel through the cassette with a hanging weight on the rack pinion emulating 15 N and 28 N excess force. Verify re-cocking at 1 mL/h cassette and the constant amplitude over the range (video amplitude tracking). Tune arm spring thickness to set E_imp.
3. **Gravimetric flow test (IEC 60601-2-24 method).** Deliver water from a 20 mL BD syringe onto a 0.001 g balance logged every 30 s for 2 h at each cassette; compute trumpet curves (2, 5, 11, 19, 31 min windows). Acceptance: mean error ≤ ±5 %, trumpet ≤ ±10 % at 2 min.
4. **Back-pressure and viscosity.** Repeat at 0, 100, 200 mmHg downstream (water column) and with 10 % dextrose at 20 °C and 37 °C; expected change < 1 % (vs 30–45 % for restriction pumps).
5. **Occlusion.** Clamp the line; record trip pressure with an in-line transducer (expect 300 ± 30 mmHg), time to trip, and bolus on release (should be < 0.2 mL).
6. **Orientation and transport.** Run on bench, tilted 90°, inverted, and on a vibrating table; rate change expected < 0.5 %.
7. **Durability.** 7-day continuous run (≈ 7 × 10⁵ band cycles) and 200 load/unload cycles of the ratchet and latch.

## 6. Limitations and honest scope

* All performance figures are **model predictions**; the models are physics-based (torsion of thin strips, beam theory, RK4 escapement dynamics, Lewis stress, Vogel viscosity) but unvalidated until step 5 is run. The quality factor of the balance in air and the printed-gear efficiencies are the least certain inputs; both were swept and both affect margin rather than rate.
* The spring-arm escapement in the CAD is a representative embodiment; its latch and arm geometry will need iteration on the bench (horological escapements always do).
* The occlusion trip must be set per syringe size (knob with marked positions); this is a training requirement.
* Rates between cassettes (e.g. 3 mL/h) require an additional cassette (any (z₁,z₂,z₃,z₄) with sums 66 can be printed).
* Pulsatility of 0.12–1.16 µL per beat is small, but very low-rate inotrope infusions (< 0.5 mL/h) would be better served with a 5 mL syringe.

## 7. Significance

ChronoFlow converts a 170-year-old horological solution to a variable-torque problem (the gravity escapement) and a precision-instrument suspension (the taut band) into a $38, printable, electricity-free infusion driver whose accuracy is limited by the syringe itself rather than by the fluid, the temperature, the back-pressure or the way it is lying. It is designed so that every safety function is mechanical and visible, and so that the rate is a physical object a nurse can hold.

---

### Files

* `chronoflow_assembly.step`, `chronoflow_assembly.glb`, `chronoflow_stl_parts.zip`, `chronoflow_viewer.html` — CAD
* `design_params.py`, `analysis_oscillator.py`, `analysis_escapement.py`, `analysis_system.py` — analysis (reproduces all figures)
* `cad_model.py`, `render_cad.py`, `check_interference.py` — model, renderer, interference check
* `results/*.csv|json` — flow table, escapement sweeps, operating map, Monte-Carlo, BOM

### Key references

AutoSyP: Juarez A. et al., *Am. J. Trop. Med. Hyg.* 2016 (doi 10.4269/ajtmh.16-0285). Electricity-free restriction pump: *ASME J. Med. Devices* 11(4):041006, 2017. Pendulum-metronome mechanical pump: ASME 2011. Spring-powered microfluidic syringe pumps: *Sens. Actuators B* 2024; *HardwareX* 2024. Gravity escapement: E. B. Denison, *Rudimentary Treatise on Clocks, Watches and Bells* (1850). Taut-band suspension: standard galvanometer practice; thin-strip torsion with axial tension: Timoshenko & Goodier, *Theory of Elasticity*. IEC 60601-2-24 (infusion pump accuracy, trumpet curves). Water viscosity: Vogel equation.

---

## Revision 2 (external review)

Four findings from external review were accepted and propagated through `design_params.py`, every analysis script, figure and document:

1. **60 mL syringe stalled before its occlusion alarm.** The first operating map was evaluated only for the 20 mL syringe. Over the full 5 syringes × 4 cassettes × 0–300 mmHg map, the 30 N spring left 3.1 N excess at the 60 mL trip (26.9 N), i.e. 3.4 µN m at the 1 mL/h cassette against a 9.5 µN m re-cock minimum — stall in 2 cells and < 1.7× in 4. The drive spring is now 2 × 5 lbf = 44.5 N; worst-cell margin 2.0× (60 mL, 1 mL/h, 300 mmHg). Consequences: torque range 19–476 µN m (24×, within the simulated spring-arm range), rack-pinion Lewis stress 22.7 MPa (SF 2.2), 4.2 J per wind, BOM $40.
2. **Trip range vs settings.** An 8–20 N knob could not deliver 6.1 or 26.9 N. The trip spring (3 N/mm) is now preloaded 2.0–9.0 mm by a spacer keyed to each syringe cradle insert. Because a force trip also senses plunger friction, the 300 mmHg setting is disclosed as a band: 280–320 mmHg (60 mL) widening to 202–398 mmHg (5 mL) for ±1.5 N friction uncertainty; the protocol measures glide force per lot.
3. **Ribbon stock.** A 0.0015" feeler gauge is 38.1 µm, not 40 µm; with k ∝ t³ that is −6.4 % in frequency, outside the trim range. The half-band is now 22 mm (was 25) for 38.1 µm, and the fixed clamp is slotted (L_h 18–26.5 mm) so that the ±2.5 µm stock tolerance (±9 % in f) is absorbed by a one-time coarse setting after measuring the ribbon. Peak band stress 166 MPa (limit 350).
4. **Mesh count.** The text said "rack 0.95, seven printed meshes at 0.92", which would give 0.95 × 0.92⁷ = 0.53. The model and CAD have seven meshes in total: the rack mesh at 0.95 and six printed gear meshes at 0.92 (four step-up stages, two cassette pairs; the escape wheel shares the cassette output shaft). 0.95 × 0.92⁶ = 0.576 is therefore the efficiency actually used throughout, and the prose was corrected to match it; no numbers changed on this account.
5. **Trip trigger travel (rev 3).** The rev-2 trip would have fired only after the full 3 mm of stem travel, adding 3 N/mm × 3 mm = 9 N to every threshold (+235 mmHg on 20 mL) and taking the 60 mL case to 9.3 µN m, below the 9.5 µN m re-cock minimum. The trip is now snap-acting: 0.3 mm sensing travel to a sear (+0.9 N, pre-compensated in the spacer preload), a separately sprung 12 mm flag rod, and sear capture of the stem so the stored compression (0.03–0.17 mL) cannot become a release bolus. Operating-map margins are unchanged.
6. **Ribbon adjustment geometry (rev 3).** Fig. 3 had shown only the upper clamp slotted. Computing k with each half's actual length, one-sided adjustment over 18–26.5 mm covers only 36.6–39.4 µm of ribbon thickness. Both clamps are now slotted and moved symmetrically (18.4–26.1 mm covers 35.6–40.6 µm).
7. **References and abstract.** Bibliography completed with full citations (Guelig 2017; Sung 2011 GHTC — previously mis-attributed to ASME IMECE; Han 2021; Park 2024 ×2; Pooke 2023; Haringx 1949; Lewis 1893; IEC 60601-2-24; ISO 7886-1). The manuscript abstract's residual "0.04 mm" was corrected to 0.0015" (38.1 µm).

## Revision 4 (external review, round 5)

Five findings accepted; all propagated through `design_params.py`, the CAD, the analyses, the figures and the documents.

8. **Occlusion mechanism was not in the CAD.** The rev-3 CAD still showed a smooth stem, a hex knob and a simple flag, so the sear, stem capture and sprung flag were claims without geometry and "zero release bolus" was an assumption. The mechanism is now fully modelled (`cad_model.trip_mechanism`, `trip_rail_linkage`): keyed preload spacer, trip spring, Ø8 stem with end collar and trip pin, sear on a vertical pivot pin and boss with return spring, Ø3 flag rod with guide boss, 2 N spring, lift and hook collars and red head, rocker (boss, cheeks, pivot pin, torsion spring) driving a Ø2 capture pin in its guide, a trip rail on two parallel-motion cranks with posts and pivot pins, and a drop rod through the deck carrying an under-deck arm and a two-pin fork around a balance spoke. The sensing travel is now derived from the drawn sear geometry (0.15 mm gap + 0.9 mm lip overlap / lever ratio 4.5 = 0.35 mm, replacing the assumed 0.3 mm; preloads 5.0–25.8 N, stored displacement 0.04–0.20 mL) and the flag travel is 6 mm (was 12). The assembly is built in armed and tripped poses, and the capture is demonstrated by solid Booleans: the stem advanced 0.1 mm after the trip interferes with the capture pin (0.085 mm³), while armed/tripped poses are clear; the fork leaves the balance free at ±0.04 rad and blocks it at ±0.08 rad, below the 6° unlock angle. Fig. 6 (schematic), Fig. 7c–d (CAD, armed vs tripped) and the text were rewritten accordingly. The arm-to-detent release link and the LOAD/RUN re-cocking of the rod and rocker remain described but not drawn.
9. **Adjustable band clamps were not in the CAD.** Both clamps were solid blocks. The upper clamp is now screwed (M3 × 2) through two slots in a 3 mm steel hanger plate that passes through a deck slot to a top flange; the lower clamp carries guide pins sliding in two slots of a bracket standing in a pocket of the bottom plate, and the tensioner's spring-seat block is screwed to the same slots, so both clamps move symmetrically over L_h = 18–26.5 mm (Fig. 7e).
10. **Constant-impulse assumption.** `analysis_arm_mechanism.py` (new Fig. 3) derives the impulse and unlatch energies from the mechanism and follows the wheel through a beat. Findings: the arm unwinds only 27 mrad per beat, so with a music-wire torsion spring (1e-4 N·m/rad, 25° preload) the impulse droops 6 % over the arc, latch-position tolerance biases it 1.4 % per build and beat-to-beat scatter is 0.07 %; a printed 0.8 mm flexure is far too soft for the preload (2.9 rad) and was replaced in the specification. Unlatch work is 0.10–0.35 µJ for µ = 0.1–0.4 (the beat model had assumed 0.12 µJ); re-running the beat simulation over E_imp ± 10 % and that range moves the period error only between 0.050 and 0.106 %. The wheel force on the locked tooth (1.6–40 mN) exceeds the arm's 8.7 mN pallet force above 105 µN·m, so a separate wheel detent is required and has been added to the CAD. The wheel-plus-reflected-cassette inertia (0.8–3.7 g·mm²) re-cocks and locks in 0.06–0.14 s against a 0.77 s budget (margin ≥ 5.5×), but the 2–98 µJ surplus per beat would load a rigid steel-pin-on-PETG detent to 2–14 N and 35–92 MPa Hertz stress (above the 50 MPa yield over most of the range, 1.7 × 10⁵ times per 20 h); the detent is therefore a 0.8 × 4 × 8 mm PETG flexure (k ≈ 2000 N/m) which limits the impact to 0.09–0.62 N and 7–19 MPa. Contact restitution and latch friction are flagged as bench inputs (protocol step 2).
11. **Interference check too permissive; non-watertight meshes.** The keyword whitelist (any pair containing "shaft", "pin", "axle", "tensioner") and the single pose were replaced by an explicit list of 20 named pairs with reasons and volume caps, evaluated in 12 poses with the gear train rotated consistently with the carriage (escape wheel snapped to a locked tooth, < 0.03 mm carriage trim). The new check found three genuine collisions that the old one could not have seen — the rack entering the cradle body at the empty-syringe position, the rail's rear post against the carriage at full travel, and the rocker against its boss in the tripped pose — all fixed (cradle rack channel, post moved to x = 238, boss lowered); it also revealed that the syringe model used the 95 mm carriage travel instead of the 69.6 mm BD 20 mL stroke (fixed) and that the flange lug was too thick for the empty-syringe thumb press (fixed). The rack polyline ran a partial last tooth past its end (self-intersecting profile); it now lays out an integer number of pitches and all 91 STLs are watertight.
12. **Reproduction.** `figstyle.py` carries verbatim copies of the plotting helpers (`apply_figure_style`, `panel_letter`, `META_GREY`, `end_of_line_labels`, `compose_figure`, `panel_px`); `analysis_system.py`, `make_panels.py`, `make_fig6.py` and `analysis_arm_mechanism.py` import them when the kernel plugins are absent, and every figure was regenerated from a plain shell. `make_renders.py` reproduces the CAD panels from the exported manifests. The README "Key numbers" are regenerated from `results/*.json` (they had been left at rev-1 values: 1.80× margin, 2.85 J, 15–275 µW, 15.3 MPa, 26.6 mm trims).
