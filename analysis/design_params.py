"""ChronoFlow -- single source of truth for the design parameters.

Electricity-free, wind-up, flexure-oscillator-governed volumetric syringe driver
for neonatal IV therapy in low-resource settings.

Units: mm, N, kg, s unless noted.
"""
import math

# ---------------------------------------------------------------- syringes (ISO 7886-1 style, BD Plastipak nominal)
SYRINGES = {
    # name: dict(ID_mm, stroke_mm (plunger travel for nominal volume), running_friction_N)
    "BD 5 mL":  dict(ID=12.06, friction=1.5),
    "BD 10 mL": dict(ID=14.50, friction=2.0),
    "BD 20 mL": dict(ID=19.13, friction=3.0),
    "BD 30 mL": dict(ID=21.70, friction=3.5),
    "BD 60 mL": dict(ID=26.70, friction=4.5),
}
for s in SYRINGES.values():
    s["area_mm2"] = math.pi * (s["ID"] / 2) ** 2
DESIGN_SYRINGE = "BD 20 mL"
SYR_AREA = SYRINGES[DESIGN_SYRINGE]["area_mm2"]  # 287.4 mm^2
SYR_VOL_ML = 20.0
SYR_STROKE = SYR_VOL_ML * 1000 / SYR_AREA  # mm of travel for 20 mL

# ---------------------------------------------------------------- energy source: constant-force (Neg'ator) spring
F_SPRING = 44.5            # N  (two 5 lbf = 22.2 N constant-force springs in parallel, e.g. Vulcan SH-series)
                           #    rev 2: raised from 30 N so that the BD 60 mL syringe at its 300 mmHg trip (26.9 N) still
                           #    leaves >= 17 N excess for the train (30 N gave 3.1 N -> escapement stalled before the trip)
F_TRIP = 14.5              # N  occlusion trip preload for the DESIGN syringe (= 300 mmHg x 287 mm2 + 3 N friction, BD 20 mL)
TRIP_SPRING_RATE = 3.0     # N/mm, trip compression spring; preload set by a keyed spacer on each syringe cradle insert
TRIP_RANGE = (6.0, 27.0)   # N  preload range the spacer set must cover (BD 5 mL ... BD 60 mL at 300 mmHg)
FRICTION_UNCERTAINTY = 1.5 # N  +/- plunger glide-force uncertainty (lot, speed, dwell) -> occlusion-pressure calibration band
def trip_force(syringe, p_mmHg=300.0):
    """Trip preload for a syringe: p x area + nominal running friction."""
    s = SYRINGES[syringe]
    return p_mmHg * 133.322 * s["area_mm2"] * 1e-6 + s["friction"]
def trip_preload(syringe, p_mmHg=300.0):
    """Spacer preload: the sear releases after TRIP_SENSE_TRAVEL, by which time the spring force has risen k*ds,
    so the preload is set k*ds BELOW the intended trip force (rev 3: the rise during sensing travel is included)."""
    return trip_force(syringe, p_mmHg) - TRIP_SPRING_RATE * TRIP_SENSE_TRAVEL
def stored_displacement_mL(syringe):
    """Pusher travel stored in the trip spring at the sear-release point, as syringe volume (upper bound on a mechanical
    release bolus if the stem were NOT captured; the sear captures the stem, so the design value is zero)."""
    return SYRINGES[syringe]["area_mm2"] * TRIP_SENSE_TRAVEL / 1000.0
def trip_pressure_range(syringe, F_set=None):
    """Pressure at which a force trip set to F_set actually fires, for friction +/- FRICTION_UNCERTAINTY (mmHg, low/nominal/high)."""
    s = SYRINGES[syringe]; F_set = trip_force(syringe) if F_set is None else F_set
    return tuple((F_set - (s["friction"] + df)) / (s["area_mm2"] * 1e-6) / 133.322 for df in (FRICTION_UNCERTAINTY, 0.0, -FRICTION_UNCERTAINTY))
# rev 4: sensing travel derived from the modelled sear geometry (cad_model.trip_mechanism):
SEAR_GAP = 0.15            # mm free gap between the stem trip pin and the sear tail (armed)
SEAR_OVERLAP = 0.9         # mm overlap of the sear hook lip over the flag-rod collar
SEAR_RATIO = 9.0 / 2.0     # hook-lip motion per tail motion (hook lever 9 mm / tail lever 2 mm about the vertical pivot)
TRIP_SENSE_TRAVEL = SEAR_GAP + SEAR_OVERLAP / SEAR_RATIO   # = 0.35 mm stem retraction at which the collar is released
TRIP_TRAVEL = TRIP_SENSE_TRAVEL
FLAG_TRAVEL = 6.0          # mm rise of the spring-loaded flag rod (cross pin enters the stem capture slot; lifts the trip rail)
STOP_FORK_FREE_PLAY_DEG = 4.8   # balance free play inside the two-pin stop fork (pins at 80/100 deg, r = 18 mm; spoke 2.5 mm wide)
SLEEVE_STOP_MM = 7.0       # lost-motion cap of the under-deck sleeve that carries the stop fork
FLAG_SPRING_N = 2.0        # N flag-rod spring (energy stored when the LOAD/RUN lever re-cocks it)
SPRING_DRUM_D = 25.0       # mm drum diameter
TRAVEL = 95.0              # mm carriage travel (covers 20 mL BD stroke 69.6 mm + loading margin)

# ---------------------------------------------------------------- drive: rack & pinion + step-up train to the escapement
MODULE = 1.0               # mm, all printed gears (FDM-friendly with 0.4 mm nozzle)
PRESSURE_ANGLE = 20.0      # deg
Z_RACK_PINION = 12         # teeth, pinion on the rack
Z_PINION = 12              # teeth, step-up pinions
Z_WHEEL = 60               # teeth, step-up wheels
STAGES = [(60, 12), (60, 12), (60, 12), (30, 12)]   # (wheel, pinion) fixed step-up stages -> 5*5*5*2.5 = 312.5
FIXED_RATIO = 1.0
for _w, _p in STAGES:
    FIXED_RATIO *= _w / _p
GEAR_FACE = 4.0            # mm face width (stages 2-4, cassette)
GEAR_FACE_STAGE1 = 8.0     # mm face width for rack, rack pinion and first wheel (highest torque)
GEAR_EFF = 0.92            # per printed mesh (lubricated PETG/PLA, lightly loaded)
RACK_EFF = 0.95

# rate cassettes: (z1, z2, z3, z4) ; speed-up ratio = (z1/z2)*(z3/z4); z1+z2 = z3+z4 = 66 (fixed centre distance 33 mm)
CASSETTES = {
    "10 mL/h": (33, 33, 33, 33),
    "5 mL/h":  (44, 22, 33, 33),
    "2 mL/h":  (33, 33, 55, 11),
    "1 mL/h":  (44, 22, 55, 11),
}
def cassette_ratio(z):
    z1, z2, z3, z4 = z
    return (z1 / z2) * (z3 / z4)

# ---------------------------------------------------------------- escapement (pin-pallet lever)
Z_ESCAPE = 15              # escape-wheel teeth
BEATS_PER_REV = 2 * Z_ESCAPE
# -- conventional pin-pallet lever escapement (baseline for comparison)
LIFT_ANGLE_DEG = 30.0      # total impulse (lift) angle seen at the balance
UNLOCK_ANGLE_DEG = 6.0     # balance angle spent unlocking (draw) before impulse
IMPULSE_FRACTION = 0.70    # fraction of the per-beat wheel energy (T_e * pi/Z) delivered as impulse
UNLOCK_FRACTION = 0.12     # fraction of per-beat wheel energy spent on unlocking (draw recoil + friction)
# -- constant-force spring-arm escapement (ChronoFlow design)
E_IMPULSE = 1.2e-6         # J per beat delivered by the cocked impulse arm (set by arm spring x impulse arc)
E_UNLATCH = 0.12e-6        # J per beat spent by the balance tripping the latch
ARM_RELEASE_DEG = 3.0      # balance angle before equilibrium at which the arm is released
ARM_IMPULSE_DEG = 15.0     # balance arc over which the arm pushes the impulse pin
LIFT_EFF = 0.6             # efficiency of the wheel re-lifting the arm (friction, impact)

# ---------------------------------------------------------------- oscillator: taut-band suspended balance
F_TARGET = None            # Hz, derived below from the flow requirement
BAND_T = 38.1e-6           # m, feeler-gauge spring steel ribbon thickness: 0.0015 in = 38.1 um (rev 2: was modelled as 40 um)
BAND_T_TOL = 2.5e-6        # m, +/- thickness tolerance of feeler stock (drives the coarse length adjustment below)
BAND_W = 2.0e-3            # m, ribbon width (cut from a 12.7 mm blade)
BAND_LH = 22.0e-3          # m, free length of each half-band (clamp to balance hub); rev 2: 25 -> 22 mm for the 38.1 um ribbon
BAND_LH_MIN, BAND_LH_MAX = 18.0e-3, 26.5e-3   # m, coarse adjustment range of BOTH end clamps (slotted posts), moved symmetrically so the
                                              #    hub stays centred; one-sided adjustment over this range would cover only 36.6-39.4 um
BAND_T0 = 1.0              # N, pretension from a steel coil spring tensioner
BAL_R_OUT = 31.0e-3        # m, balance ring outer radius
BAL_R_IN = 27.0e-3         # m
BAL_THICK = 5.0e-3         # m
RHO_PETG = 1270.0          # kg/m^3
RHO_BRASS = 8500.0
N_TRIM = 4                 # brass M3 nuts pressed into rim pockets (regulation / inertia)
M_TRIM = 0.45e-3           # kg each (M3 brass hex nut)
R_TRIM = 26.5e-3           # m  (nominal; trim nuts slide in radial slots 23-29 mm for regulation)
R_TRIM_MIN, R_TRIM_MAX = 23.0e-3, 29.0e-3
HUB_M = 1.0e-3             # kg, printed hub + impulse pin
HUB_R = 5.0e-3
THETA_MAX = 1.2            # rad, banking limit (half-amplitude)
Q_AIR = (60.0, 150.0, 400.0)   # plausible range of quality factor in air for a printed ring at 2-3 Hz

# ---------------------------------------------------------------- derived kinematics
def flow_mL_h(f_hz, cassette_ratio_value, area_mm2=SYR_AREA):
    """Volumetric flow for oscillator frequency f (Hz), cassette speed-up ratio, syringe area."""
    esc_rev_s = 2 * f_hz / BEATS_PER_REV
    pinion_rev_s = esc_rev_s / (FIXED_RATIO * cassette_ratio_value)
    v_mm_s = math.pi * MODULE * Z_RACK_PINION * pinion_rev_s
    return area_mm2 * v_mm_s * 3600 / 1000.0

def required_frequency(target_flow=10.0, cassette_ratio_value=1.0, area_mm2=SYR_AREA):
    return target_flow / flow_mL_h(1.0, cassette_ratio_value, area_mm2)

F_TARGET = required_frequency(10.0, 1.0)   # Hz -- about 1.20 Hz

OCCLUSION_LIMIT_MMHG = 300.0   # neonatal occlusion alarm limit we want the spring to enforce
MMHG_PER_PA = 1 / 133.322

if __name__ == "__main__":
    print(f"syringe area {SYR_AREA:.1f} mm2, stroke {SYR_STROKE:.1f} mm")
    print(f"fixed ratio {FIXED_RATIO:.0f}, target f {F_TARGET:.4f} Hz")
    for k, z in CASSETTES.items():
        print(k, z, cassette_ratio(z), f"{flow_mL_h(F_TARGET, cassette_ratio(z)):.3f} mL/h")
