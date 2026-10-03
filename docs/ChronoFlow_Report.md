# ChronoFlow: an electricity-free, volumetric, clockwork-governed syringe infusion driver for neonatal care in low-resource settings

**Category:** Engineering — Biomedical Engineering / Engineering Mechanics
**Status:** Complete design, parametric 3D CAD and physics-based validation. No physical prototype has been built yet; a build-and-test protocol is provided.

---

## Abstract

Hospitals without reliable electricity deliver intravenous drugs and fluids to newborns by gravity drip or by hand, with dosing errors that are often fatal at neonatal volumes. Existing non-electric infusion devices meter flow by pushing fluid through a narrow restriction, so the delivered rate depends on fluid viscosity (temperature, drug) and on downstream pressure; electric low-cost pumps still need batteries and cost hundreds of dollars. ChronoFlow is a wind-up syringe driver that meters *volume*, not pressure: a constant-force spring pushes the plunger, and a rack-and-pinion drives a 3D-printed step-up gear train into a clockwork escapement whose balance wheel sets the rate. Two design problems were identified and solved by simulation. (1) The mechanical power available to the escapement varies ~18× across clinical rates and loads, so a conventional lever escapement over-banks or stalls; adapting the 19th-century gravity-escapement principle with springs (a *constant-force spring-arm escapement*) delivers a fixed impulse energy per beat and holds the balance amplitude constant over the whole range (period error < 0.1 %). (2) Pivot friction and orientation sensitivity were eliminated with a *taut-band torsion suspension* made from a 0.04 mm feeler-gauge steel ribbon (frequency shift ±0.16 % between upright and inverted). Interchangeable colour-coded gear cassettes give 1, 2, 5 and 10 mL/h with a 20 mL syringe (0.4–19.5 mL/h across 5–60 mL syringes) at 0.12–1.16 µL per beat resolution. A 20 000-draw Monte-Carlo error budget predicts 95 % of units within ±1.14 % of set rate (requirement ±5 %), compared with +45 % at 37 °C and −33 % at 100 mmHg back-pressure for flow-restriction pumps. A calibrated trip spring stops the oscillator and raises a flag at 300 mmHg occlusion, and a syringe-latch interlock prevents free-flow boluses. The complete 66-part assembly was modelled parametrically (CadQuery), checked for interferences, and exported as STEP/STL; estimated bill of materials is $38.

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

Two constant-force springs (2 × 3.4 lbf, 30 N total) pull the carriage toward the syringe nozzle. The plunger load (syringe friction 2–5 N plus pressure × bore area) absorbs part of this force; the *excess* back-drives the rack pinion through the step-up train to the escape wheel. The escapement releases one half-tooth per beat of the balance, so plunger displacement per beat is fixed by geometry:

Δx per beat = π·m·z₀ / (2·Z_e · N_fixed · r_cassette) = π·1·12 / (2·15·312.5·1) = 4.0 µm (10 mL/h cassette), i.e. 1.16 µL per beat on a 20 mL syringe (0.12 µL at 1 mL/h).

Because the metered quantity is displacement, fluid viscosity, temperature and downstream pressure do not enter the rate (Fig. 3) until the occlusion trip point. The spring never "runs down": force is constant over the whole 95 mm travel, so runtime is set only by syringe volume (2–20 h per 20 mL fill).

### 2.2 Oscillator: taut-band balance

A 62 mm printed PETG ring (7.4 g, I = 5.2 × 10⁻⁶ kg m²) is clamped at the midpoint of a 0.04 × 2 mm spring-steel ribbon (cut from a 0.0015" feeler-gauge blade) held at both ends and pre-tensioned to 1 N by a steel coil spring. Each 25 mm half-band is a torsion spring, k_half = GJ/L + T w²/(12 L) (Timoshenko), giving k = 2.96 × 10⁻⁴ N m rad⁻¹ and f₀ = 1.20 Hz. Four brass M3 nuts in radial slots (23–29 mm) regulate the frequency by −2.2 % … +3.2 %, exactly as a clock is regulated. Design checks:

| Property | Value | Note |
|---|---|---|
| Frequency | 1.2017 Hz (target 1.2017) | trims at 26.6 mm |
| Orientation sensitivity | ±0.16 % | gravity changes tension of one half-band only |
| Peak shear stress at 1.2 rad | 154 MPa | fatigue limit ≈ 350 MPa → infinite life |
| Temperature coefficient | −1.8 × 10⁻⁴ K⁻¹ (−0.45 % over 25 K) | steel G drift + ring expansion |
| Band sag (device on its side) | 0.46 mm | within clearances |

An earlier candidate, a cross-spring flexure pivot, was rejected by analysis: a finite-element beam-column model showed that gravity on a 10 g balance changes a 0.05 mm leaf pivot's rotational stiffness by ≈20 % (the relevant scale is EI/L² ≈ 0.1 N, not the Euler load), which would make the rate orientation-dependent.

### 2.3 Escapement: constant-force spring-arm

The escapement must run from 15 µW (1 mL/h at the occlusion trip point) to 275 µW (10 mL/h, empty line) — an 18.4× torque range at the escape wheel (17–315 µN m). A balance driven by a conventional lever escapement receives an impulse proportional to that torque, so its amplitude scales as √torque: the simulation (Fig. 2a) shows it banking at > 1.2 rad over most of the operating range, which both damages rate and risks over-banking. The ChronoFlow escapement instead uses two spring-loaded impulse arms. Each escape-wheel tooth *cocks* an arm against a latch; when the balance passes, its impulse pin trips the latch and the arm delivers a fixed energy E_imp = 1.2 µJ over a 15° arc. Train torque only re-cocks the arm, so any excess is dissipated in the locking impact and never reaches the balance — the spring analogue of Denison's gravity escapement in tower clocks, whose purpose was exactly to isolate the pendulum from variable train torque.

![Figure 2 — escapement comparison]({{artifact:art_09f770dc-3d9a-4622-802c-b8ca2ef980ed}})

*Figure 2. Balance–escapement dynamics (RK4, Q = 150). (a) Amplitude vs escape-wheel torque: lever escapement amplitude grows as √T and exceeds the banking limit inside the operating range; the spring-arm escapement holds 0.59 rad from the minimum re-cock torque upward. (b) Period error relative to the free oscillator stays below 0.06 % for the spring-arm design (Airy escapement error from the asymmetric impulse window).*

Minimum re-cock torque is T_min = E_imp/(η_lift·π/Z_e) = 9.5 µN m; the worst operating point (1 mL/h at the trip point) supplies 17.1 µN m, a 1.79× margin (Fig. 5). With unknown air damping, the amplitude settles between 0.31 rad (Q = 40) and 0.99 rad (Q = 800) while the period error stays within 0.01–0.45 %, which is removed by trim calibration.

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

* **Occlusion trip.** The pusher head rides on a stem with a knob-adjustable compression spring (8–20 N). When plunger reaction exceeds the preload, the stem retracts 3 mm, exposes a red flag and lifts a trip rail that blocks the balance: the ticking stops and pressure cannot rise further. Default 14.5 N = 300 mmHg on a 20 mL syringe; settings for other syringes: 6.1 N (5 mL), 8.6 N (10 mL), 18.3 N (30 mL), 26.9 N (60 mL).
* **Free-flow protection.** The LOAD/RUN lever opens the flange latch *and* lifts the rack-pinion click at the same time. The carriage can be moved by hand only when the syringe is unlatched, so a fast-forwarded pusher moves the whole syringe along the cradle rather than expelling fluid.
* **Winding = loading.** Pulling the carriage back against the 30 N spring both loads a fresh syringe and stores 2.85 J — about eight times the 0.37 J the escapement consumes in a 20 h infusion.
* **Audible/visual running indication.** The balance beats at 1.2 Hz behind a window; silence means stopped.

### 2.6 Structural check

The highest-loaded teeth are the rack and 12-tooth rack pinion: Lewis bending stress at 30 N on an 8 mm face is 15.3 MPa against ≈50 MPa PETG yield (safety factor 3.3). Torque falls 5× per stage (180 → 33 → 6.1 → 1.1 mN m), so the remaining printed gears are lightly loaded; shafts are 3 mm steel pins in printed bores.

## 3. Predicted performance

![Figure 3 — comparison]({{artifact:art_028178a8-f664-4eff-9b94-4ec5011cb8a3}})

*Figure 3. (a) Flow-restriction pumps follow Q ∝ ΔP/μ: water viscosity (Vogel equation) falls 31 % from 20 °C to 37 °C, so a pump calibrated at 20 °C over-delivers by 45 % at body temperature and under-delivers 23 % at 10 °C; ChronoFlow drifts −0.36 % … +0.18 %. (b) A 100 mmHg venous/catheter back-pressure cuts a 300 mmHg restriction pump by 33 %; ChronoFlow is unaffected until it trips at 300 mmHg.*

![Figure 4 — error budget]({{artifact:art_cad7eac5-9c28-41b0-b52a-c9bcd2044302}})

*Figure 4. Monte-Carlo steady-state flow error (n = 20 000). Inputs: syringe bore ±1 % (uniform), trim regulation residual ±0.3 %, rack pitch after print-scale calibration ±0.2 %, ambient 15–35 °C, orientation ±0.16 %, escapement error for Q = 60–400. Result: 95 % interval −1.13 … +1.14 %, 99th percentile |error| 1.35 %. The dominant term, syringe bore tolerance, is common to every syringe pump including electronic ones.*

![Figure 5 — operating map]({{artifact:art_f7555092-cefb-404b-8274-e5a3aa1d1b45}})

*Figure 5. Escape-wheel torque relative to the arm re-cock minimum for each cassette, at the occlusion trip point (15.5 N excess) and free-running (28.5 N excess). Train efficiency 0.58 (rack 0.95 × seven printed meshes at 0.92).*

Resolution: 1.16 / 0.58 / 0.23 / 0.12 µL per beat (10/5/2/1 mL/h), i.e. 2.4 micro-boluses per second — smoother than a gravity drip (≈50 µL drops) and comparable to stepper-driven pumps.

## 4. CAD model and manufacturability

![Escapement close-up]({{artifact:art_e9ee872b-977c-49e5-878b-42357f2c15b1}})

![Exploded view]({{artifact:art_f029c561-96b4-44c1-b30f-ee6b681db7e4}})

The model is fully parametric: `design_params.py` holds every dimension and ratio; `cad_model.py` (CadQuery) builds the 66 parts in a common frame, derives shaft positions from the gear centre distances, phases every gear pair so a wheel tooth faces a pinion gap, and exports STEP, glTF and per-part STL. A pairwise solid-interference check (bounding-box prefilter + OCC Boolean common) over all 2 145 part pairs reports no unintended interferences; the only residual overlap is a 0.55 mm³ tip/root contact on the 60T/12T mesh, the known undercut of a 12-tooth 20° pinion (remedy: +0.3 profile shift). All printed parts are FDM-friendly (module 1 teeth, ≥ 2 mm walls, 0.8 mm flexures) and print in ≈420 g of PETG. Bill of materials ≈ $38 (springs $13, printed parts $9, rails $3.60, cassettes $4.80, hardware/pins/feeler stock $7.80).

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
