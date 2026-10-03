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
F_SPRING = 30.0            # N  (two 3.4 lbf constant-force springs in parallel, e.g. Vulcan SH-series)
F_TRIP = 14.5              # N  occlusion trip preload on the pusher stem (= 300 mmHg x 287 mm2 + 3 N friction, BD 20 mL); knob-adjustable 8-20 N
TRIP_TRAVEL = 3.0          # mm stem retraction that flags and stops the oscillator
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
BAND_T = 0.04e-3           # m, feeler-gauge spring steel ribbon thickness (0.0015")
BAND_W = 2.0e-3            # m, ribbon width (cut from a 12.7 mm blade)
BAND_LH = 25.0e-3          # m, free length of each half-band (clamp to balance hub)
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
