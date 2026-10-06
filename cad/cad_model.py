"""ChronoFlow parametric CAD (CadQuery).  All dimensions mm.  Global frame:
X = syringe axis / carriage travel (nozzle toward -X), Y = transverse (+Y front), Z up.
Deck plate top at Z = 0.  Upper deck: syringe cradle, carriage, rack, rack pinion.
Movement (gear train, cassette, escapement, taut-band balance) hangs below the deck
between the deck plate and a bottom plate, visible through a side window.
"""
import math, os, json
import cadquery as cq
import design_params as dp

OUT = "cad_out"; os.makedirs(OUT, exist_ok=True)
M = dp.MODULE; PA = math.radians(dp.PRESSURE_ANGLE)
SHAFT_D = 3.0; BORE = 3.15

# ------------------------------------------------------------------ layout
SYR_Y, SYR_Z = 45.0, 55.0
RAIL_Z, RAIL_Y = 46.0, (29.0, 61.0)
RACK_Y_PITCH = 18.0           # pitch line of the rack (teeth point to -Y)
RACK_Z = (26.0, 34.0)
CAR_X_MIN, CAR_X_MAX = 103.0, 213.0     # carriage centre travel (rack 110 mm always spans the pinion)
CAR_LEN = 40.0
S = {}                        # shaft plan positions
S["S1"] = (158.0, RACK_Y_PITCH - M*dp.Z_RACK_PINION/2)   # (158, 12)
def step(p, d, ang):
    return (p[0] + d*math.cos(math.radians(ang)), p[1] + d*math.sin(math.radians(ang)))
S["S2"] = step(S["S1"], 36, 200); S["S3"] = step(S["S2"], 36, -20)
S["S4"] = step(S["S3"], 36, 200); S["S5"] = step(S["S4"], 21, -20)
S["S6"] = step(S["S5"], 33, 200); S["S7"] = step(S["S6"], 33, -20)
BAL_XY = (62.0, -62.0)
# Z levels (gear mid-planes); faces
LV = dict(rack=30.0, L1=20.0, L2=12.0, L3=5.0, L4=-12.0, L5=-20.0, L6=-28.0, L7=-36.0)
DECK_Z = (-6.0, 0.0)
BOTTOM_Z = (-80.0, -74.0)
BASE_X, BASE_Y = (0.0, 260.0), (-105.0, 80.0)

COL = dict(frame="#8a9ba8", print="#2f6f9f", print2="#e8a33d", steel="#b8bcc2", brass="#c8a951",
           syringe="#d9e8f0", rubber="#333333", spring="#7f8c8d", flag="#c0392b", band="#9aa5b1")

# ------------------------------------------------------------------ gear generators
def involute_gear(z, m=M, face=4.0, bore=BORE, pa=PA, pts_per_flank=10, phase_deg=0.0):
    r_p = m*z/2; r_b = r_p*math.cos(pa); r_a = r_p + m; r_f = r_p - 1.25*m
    # involute from base circle; tooth thickness half-angle at pitch
    def inv(a): return math.tan(a) - a
    half_t = math.pi/(2*z) + inv(pa) - 0.06/r_p     # 0.06 mm tooth thinning = backlash allowance
    pts = []
    for i in range(z):
        c = 2*math.pi*i/z
        r_start = max(r_b, r_f)
        flank = []
        for j in range(pts_per_flank + 1):
            r = r_start + (r_a - r_start)*j/pts_per_flank
            a = math.acos(min(1.0, r_b/r)); phi = inv(a)
            flank.append((r, half_t - phi))
        # left flank (ccw) then right flank
        left = [(r*math.cos(c - th), r*math.sin(c - th)) for r, th in flank]
        right = [(r*math.cos(c + th), r*math.sin(c + th)) for r, th in reversed(flank)]
        # root points between teeth
        root_a = c - half_t - (math.pi/z - half_t)*0.5
        pts += [(r_f*math.cos(root_a), r_f*math.sin(root_a))] + left + right
    g = cq.Workplane("XY").polyline(pts).close().extrude(face)
    if bore: g = g.faces(">Z").workplane().hole(bore)
    if phase_deg: g = g.rotate((0, 0, 0), (0, 0, 1), phase_deg)
    return g

def mesh_phases(z_w, z_p, p_w, p_p):
    """Rotations so the wheel has a tooth centre on the line of centres and the pinion a gap there."""
    lc = math.degrees(math.atan2(p_p[1]-p_w[1], p_p[0]-p_w[0]))
    ph_w = lc % (360.0/z_w)
    ph_p = (lc + 180.0 + 180.0/z_p) % (360.0/z_p)
    return ph_w, ph_p

def rack(length, face=8.0, height=8.0, m=M, pa=PA):
    """Rack along X, teeth pointing -Y, pitch line at y=0, body from y=0 to +height.
    rev 4: an integer number of pitches is laid out symmetrically about x=0, so the outline is simple (the previous
    version ran a partial last tooth past the end and produced a self-intersecting profile -> non-watertight mesh)."""
    p = math.pi*m; ha = m; hf = 1.25*m
    t = math.tan(pa)
    n = int(length // p); L = n * p
    pts = [(-L/2, height), (-L/2, hf)]
    for i in range(n):
        xc = -L/2 + (i + 0.5) * p
        pts += [(xc - p/4 - hf*t, hf), (xc - p/4 + ha*t, -ha), (xc + p/4 - ha*t, -ha), (xc + p/4 + hf*t, hf)]
    pts += [(L/2, hf), (L/2, height)]
    return cq.Workplane("XY").polyline(pts).close().extrude(face)

def escape_wheel(z=dp.Z_ESCAPE, r=12.0, face=4.0):
    pts = []
    for i in range(z):
        a0 = 2*math.pi*i/z; a1 = a0 + 2*math.pi/z
        tip = (r*math.cos(a0), r*math.sin(a0))
        # radial locking face then sloping back to the root of the next tooth
        root = ((r-3.0)*math.cos(a0 + 0.08), (r-3.0)*math.sin(a0 + 0.08))
        back = ((r-3.0)*math.cos(a1 - 0.25), (r-3.0)*math.sin(a1 - 0.25))
        pts += [tip, root, back]
    return cq.Workplane("XY").polyline(pts).close().extrude(face).faces(">Z").workplane().hole(BORE)

def at(solid, x=0, y=0, z=0):
    return solid.translate((x, y, z))

# ------------------------------------------------------------------ frame
def deck_plate():
    p = cq.Workplane("XY").box(BASE_X[1]-BASE_X[0], BASE_Y[1]-BASE_Y[0], DECK_Z[1]-DECK_Z[0], centered=(False, False, False))
    p = p.translate((BASE_X[0], BASE_Y[0], DECK_Z[0]))
    # shaft bearing holes
    for k, (x, y) in S.items():
        p = p.cut(cq.Workplane("XY").circle(SHAFT_D/2 + 0.1).extrude(10).translate((x, y, -8)))
    # window for the balance band (upper clamp boss is separate)
    return p

def bottom_plate():
    x0, x1 = 40.0, 200.0; y0, y1 = -105.0, 30.0
    p = cq.Workplane("XY").box(x1-x0, y1-y0, BOTTOM_Z[1]-BOTTOM_Z[0], centered=(False, False, False)).translate((x0, y0, BOTTOM_Z[0]))
    for k, (x, y) in S.items():
        if k in ("S4", "S5", "S6", "S7"):
            p = p.cut(cq.Workplane("XY").circle(SHAFT_D/2 + 0.1).extrude(10).translate((x, y, BOTTOM_Z[0]-2)))
    return p

def pillars():
    ps = None
    for (x, y) in [(45, -100), (195, -100), (45, 25), (195, 25), (120, -100)]:
        c = cq.Workplane("XY").circle(4).extrude(DECK_Z[0]-BOTTOM_Z[1]).translate((x, y, BOTTOM_Z[1]))
        ps = c if ps is None else ps.union(c)
    return ps

def side_window_frame():
    # a thin acrylic window on the -Y side so the balance is visible (clinical "is it running?" cue)
    return cq.Workplane("XY").box(150, 2, 68, centered=(False, False, False)).translate((45, -107, -74))

def feet():
    f = None
    for (x, y) in [(45, -100), (195, -100), (45, 20), (195, 20)]:
        c = cq.Workplane("XY").circle(6).extrude(6).translate((x, y, BOTTOM_Z[0]-6))
        f = c if f is None else f.union(c)
    return f

# ------------------------------------------------------------------ upper deck
def rail_posts():
    posts = None
    for x in (40.0, 255.0):
        b = cq.Workplane("XY").box(8, 40, RAIL_Z + 4, centered=(True, True, False)).translate((x, SYR_Y + 1, 0))
        for y in RAIL_Y:
            b = b.cut(cq.Workplane("YZ").circle(3.1).extrude(10).translate((x-5, y, RAIL_Z)))
        b = b.cut(cq.Workplane("YZ").circle(12.5).extrude(10).translate((x-5, SYR_Y, SYR_Z)))        # barrel clearance
        b = b.cut(cq.Workplane("XY").box(10, 30, 7, centered=(True, True, False)).translate((x, SYR_Y, 29.0)))  # tape slot
        posts = b if posts is None else posts.union(b)
    return posts

def rails():
    r = None
    for y in RAIL_Y:
        c = cq.Workplane("YZ").circle(3.0).extrude(223).translate((36, y, RAIL_Z))
        r = c if r is None else r.union(c)
    return r

def cradle():
    """Syringe saddle (V-block) with flange slot at the plunger end and a tape channel underneath."""
    x0, x1 = 46.0, 98.0
    body = cq.Workplane("XY").box(x1-x0, 44, SYR_Z-8.0, centered=(False, True, False)).translate((x0, SYR_Y, 0))
    vee = cq.Workplane("YZ").polyline([(SYR_Y-24, SYR_Z+8.5), (SYR_Y, SYR_Z-15.5), (SYR_Y+24, SYR_Z+8.5)]).close().extrude(x1-x0+2).translate((x0-1, 0, 0))
    body = body.cut(vee)
    tape = cq.Workplane("XY").box(x1-x0+14, 30, 7, centered=(False, True, False)).translate((x0-1, SYR_Y, 29.0))
    body = body.cut(tape)
    for y in RAIL_Y:
        body = body.cut(cq.Workplane("YZ").circle(3.6).extrude(x1-x0+2).translate((x0-1, y, RAIL_Z)))
    # flange holder: two lugs with a 3.6 mm slot at x = 100..103.6; also carries the rails
    lug = cq.Workplane("XY").box(8.5, 36, SYR_Z+18, centered=(False, True, False)).translate((98, SYR_Y, 0))
    for y in RAIL_Y:
        lug = lug.cut(cq.Workplane("YZ").circle(3.1).extrude(12).translate((97, y, RAIL_Z)))
    slot = cq.Workplane("XY").box(3.6, 25, 40, centered=(False, True, False)).translate((100, SYR_Y, SYR_Z-18))
    lug = lug.cut(slot)
    bore = cq.Workplane("YZ").circle(11.2).extrude(12).translate((97, SYR_Y, SYR_Z))
    lug = lug.cut(bore)
    lug = lug.cut(tape)
    body = body.union(lug)
    # rev 4: clearance channel for the carriage rack (y 16.7..26, z 26..34), which reaches x = 85 with an empty syringe
    body = body.cut(cq.Workplane("XY").box(x1 - x0 + 14, 12.0, 11.0, centered=(False, False, False)).translate((x0 - 1, 15.5, 25.0)))
    return body

def barrel_clamp():
    """Hinged clamp bar over the barrel (RUN position shown)."""
    x = 60.0
    bar = cq.Workplane("XY").box(8, 52, 6, centered=(True, True, False)).translate((x, SYR_Y+3, SYR_Z+9))
    saddle = cq.Workplane("YZ").circle(11.0).extrude(10).translate((x-5, SYR_Y, SYR_Z))
    bar = bar.cut(saddle)
    hinge = cq.Workplane("XY").box(8, 6, SYR_Z+15, centered=(True, True, False)).translate((x, SYR_Y+26, 0))
    return bar.union(hinge)

def load_run_lever():
    """LOAD/RUN interlock lever on the flange holder: in LOAD it opens the flange latch AND lifts
    the rack-pinion ratchet pawl; in RUN (shown) the syringe is latched and the drive engaged."""
    base = cq.Workplane("XY").box(6, 8, 30, centered=(True, True, False)).translate((103, SYR_Y+22, SYR_Z-5))
    knob = cq.Workplane("XY").box(6, 6, 6).translate((103, SYR_Y+22, SYR_Z+27))
    return base.union(knob)

def carriage(xc):
    """Carriage body on the rails with integrated rack, pusher head, trip stem, flag."""
    body = cq.Workplane("XY").box(CAR_LEN, 50, 40, centered=(True, False, False)).translate((xc, 20, 26))
    for y in RAIL_Y:
        body = body.cut(cq.Workplane("YZ").circle(3.2).extrude(CAR_LEN+2).translate((xc-CAR_LEN/2-1, y, RAIL_Z)))
    # lighten
    body = body.cut(cq.Workplane("XY").box(CAR_LEN-8, 20, 30, centered=(True, False, False)).translate((xc, 36, 30)))
    # stem bore along X at the syringe axis
    body = body.cut(cq.Workplane("YZ").circle(4.2).extrude(CAR_LEN+2).translate((xc-CAR_LEN/2-1, SYR_Y, SYR_Z)))
    lip = cq.Workplane("XY").box(6, 28, 4, centered=(False, True, False)).translate((xc-CAR_LEN/2, SYR_Y, 33.8))
    body = body.union(lip)
    # rev 4: flag-rod and capture-pin bores through the carriage top (snap-action occlusion trip, see trip_mechanism())
    body = body.cut(cq.Workplane("XY").circle(FLAG_ROD_R + 0.1).extrude(10).translate((xc + FLAG_DX, FLAG_Y, 58.0)))
    body = body.cut(cq.Workplane("XY").circle(PIN_R + 0.1).extrude(10).translate((xc + CAP_DX, SYR_Y, 58.0)))
    rk = rack(110.0, face=RACK_Z[1]-RACK_Z[0], height=6.0)
    rk = rk.translate((xc, RACK_Y_PITCH, RACK_Z[0]))
    # rack body occupies y from pitch line to +8 -> attach to carriage underside (y 18..26)
    return body.union(rk)

# ------------------------------------------------------------------ occlusion trip (rev 4: modelled)
# Layout (all relative to the carriage centre xc; SYR_Y = 45, SYR_Z = 55, carriage top z = 66, rear face at xc+20):
#   stem        O8 steel, xc-26 .. xc+46; retracts +x against the trip spring when plunger reaction rises
#   spacer      keyed ring xc+20 .. xc+20+t_s  (t_s = preload / k sets the per-syringe trip force)
#   spring      k = 3 N/mm between spacer and end collar (xc+41 .. xc+46)
#   trip pin    O2 pin on the end collar, SEAR_GAP = 0.15 mm behind the sear tail -> sear starts turning at 0.15 mm
#   sear        2 mm plate on a vertical pivot (xc+22, 51); its hook lip overlaps the flag-rod collar by 0.9 mm;
#               lever ratio hook/tail = 9/2 = 4.5 -> lip clears the collar after a further 0.20 mm of stem travel
#               => TRIP_SENSE_TRAVEL = 0.15 + 0.9/4.5 = 0.35 mm (design_params.TRIP_SENSE_TRAVEL is derived from these)
#   flag rod    O3 rod at (xc+13, 67), 2 N spring on the guide boss; released, it rises FLAG_TRAVEL = 6 mm
#   rocker      see-saw on the carriage top (pivot y = 52): its long end is held down by the flag-rod lift collar;
#               when the rod rises, a torsion spring turns the rocker and its short end drives the capture pin down
#   capture pin O2 pin in a guide boss over the stem; drops ~3 mm into a slot in the stem top whose -x wall sits
#               exactly TRIP_SENSE_TRAVEL behind the pin, so the captured stem cannot spring forward (release bolus 0)
#   trip rail   frame-mounted bar (y = 72, z = 100) hinged at x = 250, resting on the flag head along the travel
#   drop rod    O3 rod at x = 96 pinned to the rail, through the deck; lifts a telescoping sleeve + under-deck arm
#               carrying a two-pin fork that rises around a balance spoke (lost motion: sleeve stops at the deck after 7 mm)
STEM_R = 4.0; STEM_X0, STEM_X1 = -26.0, 46.0
SPACER_R = 7.5; SPRING_R_O, SPRING_R_I = 6.0, 4.8
COLLAR_X0, COLLAR_X1, COLLAR_R = 41.0, 46.0, 6.0
SEAR_GAP = 0.15; SEAR_OVERLAP = 0.9
SEAR_PIVOT = (22.0, 51.0); SEAR_Z = (84.0, 86.0)
HOOK_TIP = (13.0, 62.0)               # hook lip centre (x rel xc, y)
TRIP_PIN_Y = 49.0                      # tail contact y -> tail lever dy = 2 mm
SEAR_RATIO = abs(HOOK_TIP[0] - SEAR_PIVOT[0]) / abs(SEAR_PIVOT[1] - TRIP_PIN_Y)   # lip y-motion per tail x-motion = 4.5
SENSE_TRAVEL = SEAR_GAP + SEAR_OVERLAP / SEAR_RATIO                               # 0.35 mm
FLAG_DX, FLAG_Y, FLAG_ROD_R = 13.0, 67.0, 1.5
FLAG_TRAVEL = 6.0
PIN_R = 1.0
CAP_DX = 6.0                           # capture pin x (rel xc) over the stem axis
ROCKER_PIVOT_Y, ROCKER_Z = 52.0, (73.5, 76.0)
ROCKER_BASE_TOP = 70.5                 # boss base top: 3 mm under the bar so the bar can tilt ~25 deg
RAIL_Y_TRIP, RAIL_ZB, RAIL_H = 72.0, 100.0, 4.0
CRANK_POSTS, CRANK_L = (70.0, 238.0), 20.0   # parallel-motion cranks: the rail translates (rises FLAG_TRAVEL everywhere)
RAIL_X0, RAIL_X1 = CRANK_POSTS[0] + CRANK_L - 2.0, CRANK_POSTS[1] + CRANK_L + 1.5
DROP_X = 96.0
STOP_R, STOP_ANGLES = 18.0, (80.0, 100.0)
DETENT_PIN_R, DETENT_ARM, DETENT_TAIL = 11.9, 12.0, 4.0   # detent pin radius on the wheel; pivot offset (tangential); L-tail length (ratio 3:1)
RESET_RAMP_X0, RESET_RAMP_X1, RESET_RAMP_GRADE = 222.0, 234.0, 1.0   # ramp underside: 6 mm drop over x = 222..228 (flag head edge reaches it at xc = 206, fully down at xc = 212), then flat to 236   # fork pins straddle the 90-degree spoke: +/-4.8 deg free play < 6 deg unlock

def _stack_rings(x0, x1, ro, ri, axis="x", yz=(0, 0), n=None):
    L = x1 - x0; n = n or max(2, int(L / 1.6))
    pitch = (L - 0.9) / max(1, n - 1)
    sp = None
    for i in range(n):
        if axis == "x":
            ring = cq.Workplane("YZ").circle(ro).circle(ri).extrude(0.9).translate((x0 + i*pitch, yz[0], yz[1]))
        else:
            ring = cq.Workplane("XY").circle(ro).circle(ri).extrude(0.9).translate((yz[0], yz[1], x0 + i*pitch))
        sp = ring if sp is None else sp.union(ring)
    return sp

def trip_mechanism(xc, tripped=False, syringe="BD 20 mL", stem_advance=0.0):
    """Snap-action occlusion trip on the carriage.  tripped=True poses the mechanism just after release:
    stem retracted by SENSE_TRAVEL, sear turned, flag rod up by FLAG_TRAVEL, rocker turned, capture pin in the stem slot."""
    d = (SENSE_TRAVEL if tripped else 0.0) - stem_advance           # stem_advance: forward creep of the stem after capture (test)
    t_s = dp.trip_preload(syringe) / dp.TRIP_SPRING_RATE            # spacer thickness (mm)
    x_face = xc - CAR_LEN/2 - 2.0
    # pusher plate + flange gripper on the stem front
    plate = cq.Workplane("YZ").rect(24, 32).extrude(4).translate((x_face-4+d, SYR_Y, SYR_Z))
    jaw = cq.Workplane("YZ").rect(24, 32).extrude(3).translate((x_face-10+d, SYR_Y, SYR_Z))
    jaw = jaw.cut(cq.Workplane("YZ").rect(18, 40).extrude(3).translate((x_face-10+d, SYR_Y, SYR_Z-10)))
    bridge = cq.Workplane("XY").box(6, 24, 4, centered=(False, True, False)).translate((x_face-10+d, SYR_Y, SYR_Z+12))
    stem = cq.Workplane("YZ").circle(STEM_R).extrude(STEM_X1 - STEM_X0).translate((xc + STEM_X0 + d, SYR_Y, SYR_Z))
    # capture slot in the stem top (armed frame, slot moves +x with the stem): the +x wall sits SENSE_TRAVEL short of the
    # pin's +x face, so after the 0.35 mm retraction it is flush with the pin and the stem cannot move back toward the
    # syringe (the pin rests on the solid stem top until then); 0.6 mm clearance on the -x side in the tripped pose
    slot_x0 = xc + CAP_DX - PIN_R - 0.6 - SENSE_TRAVEL; slot_x1 = xc + CAP_DX + PIN_R - SENSE_TRAVEL
    stem = stem.cut(cq.Workplane("XY").box(slot_x1 - slot_x0, 3.0, 4.0, centered=(False, True, False))
                    .translate((slot_x0 + d, SYR_Y, SYR_Z + STEM_R - 2.5)))
    # end collar + trip pin; spring; keyed spacer
    collar = cq.Workplane("YZ").circle(COLLAR_R).extrude(COLLAR_X1 - COLLAR_X0).translate((xc + COLLAR_X0 + d, SYR_Y, SYR_Z))
    tail_x1 = xc + COLLAR_X0 + 2.5            # sear tail +x face
    trip_pin = cq.Workplane("XY").circle(PIN_R).extrude(SEAR_Z[1] + 0.5 - (SYR_Z + 3.0)).translate(
        (tail_x1 + SEAR_GAP + PIN_R + d, TRIP_PIN_Y, SYR_Z + 3.0))
    pusher = plate.union(jaw).union(bridge).union(stem).union(collar).union(trip_pin)
    spacer = cq.Workplane("YZ").circle(SPACER_R).circle(STEM_R + 0.1).extrude(t_s).translate((xc + CAR_LEN/2, SYR_Y, SYR_Z))
    key = cq.Workplane("XY").box(t_s, 3.0, 4.0, centered=(False, True, False)).translate((xc + CAR_LEN/2, SYR_Y, SYR_Z - SPACER_R - 3.5))
    spacer = spacer.union(key)                # key tab engages a keyway in the cradle insert (one spacer per syringe size)
    spring = _stack_rings(xc + CAR_LEN/2 + t_s + 0.3, xc + COLLAR_X0 + d - 0.3, SPRING_R_O, SPRING_R_I, yz=(SYR_Y, SYR_Z))
    # ---- sear (2 mm plate on a vertical pivot pin, resting on a boss on the carriage top)
    px, py = xc + SEAR_PIVOT[0], SEAR_PIVOT[1]
    hub = cq.Workplane("XY").circle(3.0).extrude(SEAR_Z[1] - SEAR_Z[0]).translate((px, py, SEAR_Z[0]))
    hx, hy = xc + HOOK_TIP[0], HOOK_TIP[1]
    L_hook = math.hypot(hx - px, hy - py); a_hook = math.degrees(math.atan2(hy - py, hx - px))
    hook_arm = (cq.Workplane("XY").box(L_hook, 3.0, SEAR_Z[1] - SEAR_Z[0], centered=(False, True, False))
                .rotate((0, 0, 0), (0, 0, 1), a_hook).translate((px, py, SEAR_Z[0])))
    # hook lip: pad whose +y edge overlaps the collar (collar -y edge at FLAG_Y-3) by SEAR_OVERLAP
    lip = cq.Workplane("XY").box(3.0, 2.4, SEAR_Z[1] - SEAR_Z[0], centered=(True, False, False)).translate(
        (hx, FLAG_Y - 3.0 + SEAR_OVERLAP - 2.4, SEAR_Z[0]))
    tail = cq.Workplane("XY").box(tail_x1 - px, 3.0, SEAR_Z[1] - SEAR_Z[0], centered=(False, True, False)).translate(
        (px, TRIP_PIN_Y, SEAR_Z[0]))
    sear = hub.union(hook_arm).union(lip).union(tail)
    sear = sear.cut(cq.Workplane("XY").circle(0.8).extrude(5).translate((px, py, SEAR_Z[0] - 1)))
    if tripped:
        psi = math.degrees((SENSE_TRAVEL - SEAR_GAP) / (py - TRIP_PIN_Y))   # tail pushed +x by (travel - gap); tail lever = py - TRIP_PIN_Y
        sear = sear.rotate((px, py, 0), (px, py, 1), psi)
    sear_pin = cq.Workplane("XY").circle(0.75).extrude(SEAR_Z[1] + 1 - 66.0).translate((px, py, 66.0))
    sear_boss = cq.Workplane("XY").box(6, 6, SEAR_Z[0] - 66.0, centered=(True, True, False)).translate((px, py, 66.0))
    sear_boss = sear_boss.cut(cq.Workplane("XY").circle(0.75).extrude(30).translate((px, py, 66.0)))
    sear_spring = cq.Workplane("XY").circle(2.2).circle(1.6).extrude(1.5).translate((px, py, SEAR_Z[1]))   # torsion-spring stand-in
    # ---- flag rod: guide boss on the carriage top, spring, lift collar (rocker), hook collar (sear), head
    dz = FLAG_TRAVEL if tripped else 0.0
    fx = xc + FLAG_DX
    rod = cq.Workplane("XY").circle(FLAG_ROD_R).extrude(100.0 - 60.0).translate((fx, FLAG_Y, 60.0 + dz))
    lift_collar = cq.Workplane("XY").circle(3.5).extrude(2.0).translate((fx, FLAG_Y, ROCKER_Z[1] + dz))          # 76..78
    hook_collar = cq.Workplane("XY").circle(3.0).extrude(2.0).translate((fx, FLAG_Y, SEAR_Z[0] - 2.0 + dz))      # 80..82
    head = cq.Workplane("XY").box(6, 16, 4, centered=(True, True, False)).translate((fx, FLAG_Y, 96.0 + dz))   # 6 mm in x: its leading edge rides the reset ramp
    flag_rod = rod.union(lift_collar).union(hook_collar).union(head)
    boss = cq.Workplane("XY").circle(4.0).circle(FLAG_ROD_R + 0.1).extrude(4.0).translate((fx, FLAG_Y, 66.0))
    flag_spring = _stack_rings(70.3, ROCKER_Z[1] + dz - 0.3, 2.8, 2.0, axis="z", yz=(fx, FLAG_Y), n=3)
    # ---- rocker (bar in the YZ plane at x = xc+CAP_DX): long end under the lift collar, short end under the capture-pin head
    rx = xc + CAP_DX
    y0, y1 = SYR_Y - 2.0, FLAG_Y - 2.0
    bar = cq.Workplane("XY").box(4.4, y1 - y0, ROCKER_Z[1] - ROCKER_Z[0], centered=(True, False, False)).translate((rx, y0, ROCKER_Z[0]))
    zp = (ROCKER_Z[0] + ROCKER_Z[1]) / 2
    bar = bar.cut(cq.Workplane("YZ").circle(0.8).extrude(6).translate((rx - 3, ROCKER_PIVOT_Y, zp)))
    bar = bar.cut(cq.Workplane("XY").box(2*PIN_R + 0.4, SYR_Y + 2.6 - (y0 - 1), 6, centered=(True, False, False))
                  .translate((rx, y0 - 1, ROCKER_Z[0] - 1)))     # fork around the capture-pin shaft (bar end sits under its head)
    cap_len = ROCKER_PIVOT_Y - SYR_Y; lift_len = FLAG_Y - ROCKER_PIVOT_Y
    phi = math.degrees(math.atan2(dz, lift_len)) if tripped else 0.0
    # pin head rests on the highest point of the tilted bar under it (bar edge at y = SYR_Y + 2): ~2.2 mm travel for a 6 mm rod rise
    cap_drop = ROCKER_Z[1] - (zp + (SYR_Y + 2.0 - ROCKER_PIVOT_Y) * math.sin(math.radians(phi)) + (ROCKER_Z[1] - zp) * math.cos(math.radians(phi)))
    if tripped:
        bar = bar.rotate((rx - 5, ROCKER_PIVOT_Y, zp), (rx + 5, ROCKER_PIVOT_Y, zp), phi)   # +phi lifts the +y (collar) end
    rocker_boss = cq.Workplane("XY").box(8, 6, ROCKER_BASE_TOP - 66.0, centered=(True, True, False)).translate((rx, ROCKER_PIVOT_Y, 66.0))
    for sx in (-2.95, 2.95):   # cheeks carrying the pivot pin
        rocker_boss = rocker_boss.union(cq.Workplane("XY").box(1.5, 6, ROCKER_Z[1] + 1.0 - ROCKER_BASE_TOP, centered=(True, True, False))
                                        .translate((rx + sx, ROCKER_PIVOT_Y, ROCKER_BASE_TOP)))
    rocker_boss = rocker_boss.cut(cq.Workplane("YZ").circle(0.8).extrude(8).translate((rx - 4, ROCKER_PIVOT_Y, zp)))
    rocker_pin = cq.Workplane("YZ").circle(0.75).extrude(8).translate((rx - 4, ROCKER_PIVOT_Y, zp))
    rocker_spring = cq.Workplane("YZ").circle(2.4).circle(1.7).extrude(1.5).translate((rx + 3.9, ROCKER_PIVOT_Y, zp))  # torsion-spring stand-in
    # ---- capture pin in its guide boss (tip 1 mm above the stem top when armed)
    cap_tip = SYR_Z + STEM_R + 0.5 - cap_drop                        # 0.5 mm above the stem top when armed
    cap_pin = cq.Workplane("XY").circle(PIN_R).extrude(ROCKER_Z[1] - cap_tip).translate((rx, SYR_Y, cap_tip))
    cap_pin = cap_pin.union(cq.Workplane("XY").circle(2.0).extrude(2.0).translate((rx, SYR_Y, ROCKER_Z[1] - cap_drop)))  # head
    cap_boss = cq.Workplane("XY").circle(3.0).circle(PIN_R + 0.1).extrude(ROCKER_BASE_TOP - 0.5 - 66.0).translate((rx, SYR_Y, 66.0))
    return dict(pusher=pusher, spacer=spacer, trip_spring=spring, sear=sear, sear_pin=sear_pin, sear_boss=sear_boss,
                sear_spring=sear_spring, flag_rod=flag_rod, flag_boss=boss, flag_spring=flag_spring,
                rocker=bar, rocker_boss=rocker_boss, rocker_pin=rocker_pin, rocker_spring=rocker_spring,
                capture_pin=cap_pin, capture_boss=cap_boss), t_s

def trip_rail_linkage(xc, tripped=False):
    """Frame-mounted trip rail on two parallel-motion cranks + drop rod + under-deck arm with a two-pin fork that rises
    around a balance spoke.  The rail is a parallelogram coupler, so it rises by FLAG_TRAVEL wherever the flag lifts it."""
    bx, by = BAL_XY
    zp = RAIL_ZB + RAIL_H/2
    alpha = math.degrees(math.asin(FLAG_TRAVEL / CRANK_L)) if tripped else 0.0
    dx = CRANK_L * (math.cos(math.radians(alpha)) - 1.0); dz = FLAG_TRAVEL if tripped else 0.0
    posts = None; cranks = None; pins = None
    for xp in CRANK_POSTS:
        post = cq.Workplane("XY").box(6, 8, RAIL_ZB - 6.0, centered=(True, True, False)).translate((xp, RAIL_Y_TRIP, 0))
        cheek = cq.Workplane("XY").box(6, 2, 10, centered=(True, True, False)).translate((xp, RAIL_Y_TRIP - 6, RAIL_ZB - 6.0))
        cheek = cheek.cut(cq.Workplane("XZ").circle(1.1).extrude(4).translate((xp, RAIL_Y_TRIP - 4, zp)))
        post = post.union(cheek)
        posts = post if posts is None else posts.union(post)
        crank = cq.Workplane("XY").box(CRANK_L + 6, 2, 4, centered=(False, True, False)).translate((xp - 3, RAIL_Y_TRIP - 4, zp - 2))
        crank = crank.cut(cq.Workplane("XZ").circle(1.1).extrude(4).translate((xp, RAIL_Y_TRIP - 2, zp)))
        crank = crank.cut(cq.Workplane("XZ").circle(1.1).extrude(4).translate((xp + CRANK_L, RAIL_Y_TRIP - 2, zp)))
        if tripped: crank = crank.rotate((xp, 0, zp), (xp, 1, zp), alpha)
        cranks = crank if cranks is None else cranks.union(crank)
        p1 = cq.Workplane("XZ").circle(1.0).extrude(4.5).translate((xp, RAIL_Y_TRIP - 3.0, zp))                      # post pivot (cheek + crank)
        p2 = cq.Workplane("XZ").circle(1.0).extrude(9).translate((xp + CRANK_L + dx, RAIL_Y_TRIP + 3.5, zp + dz))    # rail pivot (crank + rail)
        pins = p1.union(p2) if pins is None else pins.union(p1).union(p2)
    rail = cq.Workplane("XY").box(RAIL_X1 - RAIL_X0, 6, RAIL_H, centered=(False, True, False)).translate((RAIL_X0, RAIL_Y_TRIP, RAIL_ZB))
    for xp in CRANK_POSTS:
        rail = rail.cut(cq.Workplane("XZ").circle(1.1).extrude(12).translate((xp + CRANK_L, RAIL_Y_TRIP + 6, zp)))
    rail = rail.cut(cq.Workplane("XY").circle(FLAG_ROD_R + 0.1).extrude(RAIL_H + 2).translate((DROP_X, RAIL_Y_TRIP, RAIL_ZB - 1)))   # drop rod seat
    rail = rail.translate((dx, 0, dz))
    # drop rod fixed to the rail, through the deck, carrying the under-deck arm and the stop fork
    y_fork = by + STOP_R * math.sin(math.radians(STOP_ANGLES[0]))
    z_arm = -50.0
    rod = cq.Workplane("XY").circle(FLAG_ROD_R).extrude(RAIL_ZB + RAIL_H - z_arm).translate((DROP_X, RAIL_Y_TRIP, z_arm))
    arm_y = cq.Workplane("XY").box(4, RAIL_Y_TRIP - y_fork + 4, 4, centered=(True, False, False)).translate((DROP_X, y_fork - 2, z_arm))
    xs = [bx + STOP_R * math.cos(math.radians(a)) for a in STOP_ANGLES]
    arm_x = cq.Workplane("XY").box(DROP_X + 2 - (min(xs) - 2), 4, 4, centered=(False, True, False)).translate((min(xs) - 2, y_fork, z_arm))
    link = rod.union(arm_y).union(arm_x)
    dx_trip = CRANK_L * (math.cos(math.asin(FLAG_TRAVEL / CRANK_L)) - 1.0)   # parallelogram side-shift when tripped (-0.93 mm)
    for a in STOP_ANGLES:   # fork pins: tips at z = -41 armed (3.75 mm under the spoke), -35 tripped (inside the 2.5 mm spoke);
        # placed so that they are symmetric about the spoke in the TRIPPED pose (pre-offset by -dx_trip)
        p = cq.Workplane("XY").circle(PIN_R).extrude(-41.0 - (z_arm + 4)).translate(
            (bx + STOP_R * math.cos(math.radians(a)) - dx_trip, by + STOP_R * math.sin(math.radians(a)), z_arm + 4))
        link = link.union(p)
    link = link.translate((dx, 0, dz))
    # rev 5: RESET RAMP.  Pulling the carriage fully back to load (xc -> CAR_X_MAX) runs the flag head under a frame-mounted
    # ramp that presses the rod down 6 mm: the hook collar re-engages the sear (sear snaps in on its return spring, collar
    # bottom chamfered) and the lift collar turns the rocker back, raising the capture pin out of the stem slot.
    ramp_x0, ramp_x1, ramp_slope = RESET_RAMP_X0, RESET_RAMP_X0 + FLAG_TRAVEL / RESET_RAMP_GRADE, RESET_RAMP_X1
    z_hi, z_lo = RAIL_ZB + FLAG_TRAVEL, RAIL_ZB                    # head top tripped -> head top armed
    prof = [(ramp_x0, z_hi), (ramp_x1, z_lo), (ramp_slope, z_lo), (ramp_slope, z_lo + 3.0), (ramp_x1, z_lo + 3.0), (ramp_x0, z_hi + 3.0)]
    ramp = cq.Workplane("XZ").polyline(prof).close().extrude(-6.0).translate((0, FLAG_Y - 6.0, 0))      # y 61..67 (beside the rail, y 69..75)
    web = cq.Workplane("XY").box(3, 6, 111.0 - (z_lo + 3.0), centered=(False, True, False)).translate((ramp_slope - 7.0, FLAG_Y - 3.0, z_lo + 3.0))
    xp = CRANK_POSTS[1]
    bracket = cq.Workplane("XY").box(xp + 3 - (ramp_slope - 7.0), 80.0 - (FLAG_Y - 6.0), 3, centered=(False, False, False)).translate((ramp_slope - 7.0, FLAG_Y - 6.0, 111.0))
    column = cq.Workplane("XY").box(6, 3.5, 111.0, centered=(True, False, False)).translate((xp, 76.5, 0))       # on the deck, outboard of the rail
    ramp_post = web.union(bracket).union(column)
    return dict(trip_rail=rail, rail_posts=posts, rail_cranks=cranks, rail_pins=pins, stop_linkage=link, reset_ramp=ramp.union(ramp_post))

def spring_drum():
    cx, cz = 18.0, 21.0
    drum = cq.Workplane("XZ").circle(dp.SPRING_DRUM_D/2).circle(2.1).extrude(28).translate((cx, SYR_Y+14, cz))
    coil = cq.Workplane("XZ").circle(dp.SPRING_DRUM_D/2 + 4).extrude(26).translate((cx, SYR_Y+13, cz))
    coil = coil.cut(cq.Workplane("XZ").circle(dp.SPRING_DRUM_D/2 + 0.05).extrude(30).translate((cx, SYR_Y+15, cz)))
    bracket = None
    for y in (SYR_Y - 18.0, SYR_Y + 15.0):
        plate = cq.Workplane("XY").box(16, 3, cz + 6, centered=(True, False, False)).translate((cx, y, 0))
        plate = plate.cut(cq.Workplane("XZ").circle(2.1).extrude(10).translate((cx, y + 6, cz)))
        bracket = plate if bracket is None else bracket.union(plate)
    bracket = bracket.union(cq.Workplane("XY").box(16, 36, 3, centered=(True, True, False)).translate((cx, SYR_Y, 0)))
    axle = cq.Workplane("XZ").circle(2.0).extrude(44).translate((cx, SYR_Y+22, cz))
    return drum, coil, bracket, axle

def spring_tape(xc):
    x0 = 18.0 + dp.SPRING_DRUM_D/2 + 0.5
    x1 = xc - CAR_LEN/2
    return cq.Workplane("XY").box(x1-x0, 26, 0.3, centered=(False, True, False)).translate((x0, SYR_Y, 33.5))

def syringe(fill_fraction=1.0):
    """Simplified BD 20 mL: barrel, flange at x=100..103, plunger, thumb press."""
    ID = dp.SYRINGES["BD 20 mL"]["ID"]; OD = ID + 2.4
    x_fl = 100.0
    barrel = cq.Workplane("YZ").circle(OD/2).circle(ID/2).extrude(88).translate((x_fl-88, SYR_Y, SYR_Z))
    end = cq.Workplane("YZ").circle(OD/2).extrude(3).translate((x_fl-91, SYR_Y, SYR_Z))
    nozzle = cq.Workplane("YZ").circle(2.2).extrude(10).translate((x_fl-101, SYR_Y, SYR_Z))
    flange = cq.Workplane("YZ").rect(24, 33).extrude(3).translate((x_fl, SYR_Y, SYR_Z))
    flange = flange.cut(cq.Workplane("YZ").circle(ID/2).extrude(5).translate((x_fl-1, SYR_Y, SYR_Z)))
    fluid_len = (dp.SYR_VOL_ML * 1000.0 / dp.SYRINGES["BD 20 mL"]["area_mm2"]) * fill_fraction   # 69.6 mm stroke of the BD 20 mL (rev 4; was the 95 mm carriage travel)
    stopper_x = x_fl - 88 + 3 + fluid_len               # front face of the rubber stopper
    rod_len = 86.0                                      # empty syringe -> thumb press 6 mm behind the flange (clears the 6.5 mm flange lug)
    stopper = cq.Workplane("YZ").circle(ID/2 - 0.1).extrude(8).translate((stopper_x, SYR_Y, SYR_Z))
    rod = cq.Workplane("YZ").rect(9, 9).extrude(rod_len).translate((stopper_x + 8, SYR_Y, SYR_Z))
    thumb_x = stopper_x + 8 + rod_len
    thumb = cq.Workplane("YZ").circle(14).extrude(3).translate((thumb_x, SYR_Y, SYR_Z))
    return barrel.union(end).union(nozzle).union(flange), stopper, rod.union(thumb), thumb_x

# ------------------------------------------------------------------ movement
def shafts():
    sh = {}
    zr = dict(S1=(10, 36), S2=(4, 26), S3=(-16, 16), S4=(-24, 10), S5=(-26, -10), S6=(-34, -16), S7=(-42, -16))
    for k, (x, y) in S.items():
        z0, z1 = zr[k]
        sh[k] = cq.Workplane("XY").circle(SHAFT_D/2).extrude(z1-z0).translate((x, y, z0))
    return sh

XC_REF = 213.0 - math.pi/2 - 17*math.pi      # carriage x at which a rack tooth is centred under the pinion (= 158.02)

def train_angles(cassette, xc):
    """Rotation (deg) of every gear for carriage position xc, consistent through all seven meshes; the carriage
    displacement is trimmed by < 0.03 mm so that the escape wheel sits at a locked tooth (multiple of 24 deg)."""
    z1, z2, z3, z4 = dp.CASSETTES[cassette]
    r_p0 = M * dp.Z_RACK_PINION / 2
    ratios = [1.0, -60/12, -60/12, -60/12, -30/12, -z1/z2, -z3/z4]       # P0/W1, P1/W2, P2/W3, P3/W4, P4/g1, g2/g3, g4/esc
    K = -1.0 / r_p0                                                    # rad of P0 per mm of carriage (+x rack -> cw)
    for r in ratios[1:]: K *= r
    dx = xc - XC_REF
    esc = math.degrees(K * dx); esc_snap = round(esc / 24.0) * 24.0
    dx = dx + (math.radians(esc_snap) - math.radians(esc)) / K if K else dx
    th = -dx / r_p0; out = []
    for r in ratios:
        th = th * r if r != 1.0 else th
        out.append(math.degrees(th))
    return dict(P0=out[0], W1=out[0], P1=out[1], W2=out[1], P2=out[2], W3=out[2], P3=out[3], W4=out[3],
                P4=out[4], g1=out[4], g2=out[5], g3=out[5], g4=out[6], esc=out[6])

def gear_train(cassette, xc=None):
    z1, z2, z3, z4 = dp.CASSETTES[cassette]
    G = {}
    f8 = dp.GEAR_FACE_STAGE1; f4 = dp.GEAR_FACE
    ang = train_angles(cassette, xc) if xc is not None else {k: 0.0 for k in "P0 W1 P1 W2 P2 W3 P3 W4 P4 g1 g2 g3 g4 esc".split()}
    # rack pinion: at XC_REF a rack tooth is centred under the pinion and a pinion gap faces +Y (the rack)
    G["P0 rack pinion"] = at(involute_gear(dp.Z_RACK_PINION, face=f8, phase_deg=(90 + 180/12) % 30 + ang["P0"]), *S["S1"], RACK_Z[0])
    pairs = [("W1 60T", 60, "S1", "P1 12T", 12, "S2", f8, LV["L1"], "W1", "P1"), ("W2 60T", 60, "S2", "P2 12T", 12, "S3", f4, LV["L2"], "W2", "P2"),
             ("W3 60T", 60, "S3", "P3 12T", 12, "S4", f4, LV["L3"], "W3", "P3"), ("W4 30T", 30, "S4", "P4 12T", 12, "S5", f4, LV["L4"], "W4", "P4"),
             (f"cassette g1 {z1}T", z1, "S5", f"cassette g2 {z2}T", z2, "S6", f4, LV["L5"], "g1", "g2"),
             (f"cassette g3 {z3}T", z3, "S6", f"cassette g4 {z4}T", z4, "S7", f4, LV["L6"], "g3", "g4")]
    for nw, zw, sw, npn, zp, sp, f, lv, kw, kp in pairs:
        phw, php = mesh_phases(zw, zp, S[sw], S[sp])
        G[nw] = at(involute_gear(zw, face=f, phase_deg=phw + ang[kw]), *S[sw], lv - f/2)
        G[npn] = at(involute_gear(zp, face=f, phase_deg=php + ang[kp]), *S[sp], lv - f/2)
    G["escape wheel 15T"] = at(escape_wheel(face=f4).rotate((0, 0, 0), (0, 0, 1), ang["esc"] % 24.0), *S["S7"], LV["L7"]-f4/2)
    return G

def ratchet_click():
    """One-way click on S1 between P0 and W1: retracting the carriage (winding) freewheels."""
    x, y = S["S1"]
    wheel = cq.Workplane("XY").polygon(12, 16).extrude(1.6).faces(">Z").workplane().hole(BORE).translate((x, y, RACK_Z[0]-1.8))
    pawl = cq.Workplane("XY").box(10, 3, 1.6, centered=(False, True, False)).translate((x+8, y+2, RACK_Z[0]-1.8))
    return wheel.union(pawl)

def cassette_plate():
    """Carrier plate that holds the idler pin S6 and keys g1/g4 to shafts S5/S7; colour-coded per rate."""
    (x5, y5), (x6, y6), (x7, y7) = S["S5"], S["S6"], S["S7"]
    z = LV["L6"]-dp.GEAR_FACE/2-3.5
    pl = None
    for (x, y) in (S["S5"], S["S6"], S["S7"]):
        d = cq.Workplane("XY").circle(7).extrude(3).translate((x, y, z))
        pl = d if pl is None else pl.union(d)
    for (a, b) in ((S["S5"], S["S6"]), (S["S6"], S["S7"])):
        L = math.hypot(b[0]-a[0], b[1]-a[1]); ang = math.degrees(math.atan2(b[1]-a[1], b[0]-a[0]))
        bar = cq.Workplane("XY").box(L, 8, 3, centered=(False, True, False)).rotate((0,0,0),(0,0,1), ang).translate((a[0], a[1], z))
        pl = pl.union(bar)
    for (x, y) in (S["S5"], S["S6"], S["S7"]):
        pl = pl.cut(cq.Workplane("XY").circle(SHAFT_D/2+0.1).extrude(5).translate((x, y, LV["L6"]-dp.GEAR_FACE/2-4)))
    return pl

def spring_arm_escapement():
    """Constant-force spring-arm escapement: two impulse arms pivoted beside the escape wheel,
    each with a steel pallet pin that is lifted (cocked) by the wheel teeth, a printed flexure
    spring, and a latch released by the balance impulse pin.  Representative embodiment."""
    ex, ey = S["S7"]; bx, by = BAL_XY
    parts = {}; RELEASE_TAILS = []
    # arm pivots: between wheel and balance, 14 mm from wheel centre toward the balance, +/-9 mm sideways
    ux, uy = (bx-ex), (by-ey); L = math.hypot(ux, uy); ux, uy = ux/L, uy/L
    nx, ny = -uy, ux
    for sgn, name, z0 in ((+1, "impulse arm A", LV["L7"] - 8.5), (-1, "impulse arm B", LV["L7"] - 5.5)):
        px, py = ex + 14*ux + sgn*9*nx, ey + 14*uy + sgn*9*ny
        # arm from pivot toward the balance axis (impulse tip near the balance hub)
        tip_x, tip_y = bx - 9*ux, by - 9*uy
        arm = cq.Workplane("XY").polyline([(px, py), (tip_x + sgn*3.0*nx, tip_y + sgn*3.0*ny), (tip_x + sgn*0.5*nx, tip_y + sgn*0.5*ny)]).close().extrude(2.5).translate((0, 0, z0))
        arm = arm.union(cq.Workplane("XY").circle(3).extrude(2.5).translate((px, py, z0)))
        # pallet pin toward the wheel rim
        wx, wy = ex + 11.5*math.cos(math.atan2(py-ey, px-ex)), ey + 11.5*math.sin(math.atan2(py-ey, px-ex))
        arm = arm.union(cq.Workplane("XY").polyline([(px, py), (wx + 1*nx*sgn, wy + 1*ny*sgn), (wx - 1*nx*sgn, wy - 1*ny*sgn)]).close().extrude(2.5).translate((0, 0, z0)))
        parts[name] = arm
        parts[name + " pallet pin"] = cq.Workplane("XY").circle(0.5).extrude(LV["L7"] + 2.0 - (z0 - 0.5)).translate((wx, wy, z0-0.5))
        parts[name + " pivot pin"] = cq.Workplane("XY").circle(SHAFT_D/2 - 0.5).extrude(12).translate((px, py, z0-7))
        # rev 5: arm preload spring = music-wire torsion spring around the pivot pin (k ~ 1e-4 N m/rad, 25 deg preload);
        # the printed 0.8 mm flexure of earlier revisions (k = 0.0149 N m/rad) is 150x too stiff to act as a constant-force
        # element (analysis_arm_mechanism.py).  Stand-in: coil above the arm, straight leg anchored on the frame post.
        sx, sy = px + sgn*8*nx, py + sgn*8*ny
        parts[name + " torsion spring"] = (_stack_rings(z0 + 2.6, z0 + 6.6, 2.6, 1.9, axis="z", yz=(px, py), n=4)
            .union(cq.Workplane("XY").box(math.hypot(sx - px, sy - py), 0.8, 0.8, centered=(False, True, False))
                   .rotate((0, 0, 0), (0, 0, 1), math.degrees(math.atan2(sy - py, sx - px))).translate((px, py, z0 + 6.6))))
        parts[name + " post"] = cq.Workplane("XY").circle(2).extrude(z0 + 8 - BOTTOM_Z[1]).translate((sx, sy, BOTTOM_Z[1]))
        # rev 5: release tail behind the pivot.  Tail A ends 1 mm beyond the wheel axis (16 mm), tail B 1 mm short of it
        # (14 mm), so that the two tails -- which move in opposite n-directions -- turn the release star the SAME way.
        Lt = 15.0 - sgn * 1.0
        tx, ty = px - Lt*ux, py - Lt*uy
        parts[name] = parts[name].union(cq.Workplane("XY").box(Lt, 2.4, 2.5, centered=(False, True, False))
                                        .rotate((0, 0, 0), (0, 0, 1), math.degrees(math.atan2(-uy, -ux))).translate((px, py, z0)))
        RELEASE_TAILS.append((sgn, tx, ty))
    # latch block between the arms, released by the balance impulse pin
    lx, ly = ex + 30*ux, ey + 30*uy
    parts["latch"] = cq.Workplane("XY").box(4, 10, 5.5, centered=(True, True, False)).rotate((0,0,0),(0,0,1), math.degrees(math.atan2(uy, ux))).translate((lx, ly, LV["L7"] - 8.5))
    # rev 4/5: compliant wheel DETENT.  The wheel is locked by this detent (not by the arm) during impulse.  rev 5 makes the
    # release path explicit: each impulse arm carries a tail that, at the end of the impulse stroke, pushes a peg of a
    # RELEASE STAR turning freely on a stub axle coaxial with the escape shaft (below the wheel); the star's third arm
    # pushes the L-tail of the detent lever, whose pivot sits tangentially beside the pin so the pin withdraws RADIALLY
    # out of the tooth gap.  The detent flexure (0.8 x 4 x 8 mm PETG, k ~ 2000 N/m) returns it and takes the lock impact.
    zd = LV["L7"] - 5.5                                               # lever -41.5..-39 (below the wheel plane); only the pin enters the teeth
    ang_back = math.atan2(-uy, -ux); pitch = 2*math.pi/dp.Z_ESCAPE
    a_pin = round(ang_back / pitch) * pitch + 0.21                    # in the gap behind a tooth root
    r_pin = DETENT_PIN_R
    rhx, rhy = math.cos(a_pin), math.sin(a_pin); thx, thy = -math.sin(a_pin), math.cos(a_pin)   # radial / tangential unit vectors at the pin
    pxd, pyd = ex + r_pin*rhx, ey + r_pin*rhy
    hub_x, hub_y = pxd + DETENT_ARM*thx, pyd + DETENT_ARM*thy          # pivot tangentially beside the pin -> radial pin motion
    tlx, tly = hub_x - DETENT_TAIL*rhx, hub_y - DETENT_TAIL*rhy         # L-tail end (toward the wheel axis)
    def bar(x0, y0, x1, y1, w, h, z):
        L_ = math.hypot(x1 - x0, y1 - y0)
        return cq.Workplane("XY").box(L_, w, h, centered=(False, True, False)).rotate((0, 0, 0), (0, 0, 1), math.degrees(math.atan2(y1 - y0, x1 - x0))).translate((x0, y0, z))
    lever = bar(hub_x, hub_y, pxd, pyd, 2.4, 2.5, zd).union(bar(hub_x, hub_y, tlx, tly, 2.4, 2.5, zd))
    lever = lever.union(cq.Workplane("XY").circle(2.5).extrude(2.5).translate((hub_x, hub_y, zd)))
    lever = lever.cut(cq.Workplane("XY").circle(1.05).extrude(4).translate((hub_x, hub_y, zd - 0.5)))
    parts["wheel detent"] = lever
    parts["wheel detent pin"] = cq.Workplane("XY").circle(0.5).extrude(LV["L7"] + 2.0 - (zd - 0.5)).translate((pxd, pyd, zd - 0.5))
    parts["wheel detent pivot pin"] = cq.Workplane("XY").circle(1.0).extrude(zd + 4 - BOTTOM_Z[1]).translate((hub_x, hub_y, BOTTOM_Z[1]))
    fx_, fy_ = hub_x + 5*rhx + 4*thx, hub_y + 5*rhy + 4*thy             # flexure outboard of the hub, pressing the lever toward the wheel
    parts["wheel detent flexure"] = cq.Workplane("XY").box(8, 0.8, 4, centered=(True, True, False)).rotate((0,0,0),(0,0,1), math.degrees(math.atan2(rhy, rhx))).translate((fx_, fy_, zd - 0.5))
    parts["wheel detent post"] = cq.Workplane("XY").circle(2).extrude(zd + 4 - BOTTOM_Z[1]).translate((fx_ + 4*rhx, fy_ + 4*rhy, BOTTOM_Z[1]))
    # release star on a stub axle coaxial with S7, below the shaft end
    zs = LV["L7"] - 11.5                                              # star plate -47.5..-45.5
    star = cq.Workplane("XY").circle(3.5).extrude(2.0).translate((ex, ey, zs))
    pegs = None
    for sgn, tx, ty in RELEASE_TAILS:                                  # pegs beside each tail end (tail pushes the peg sideways)
        qx, qy = tx + sgn*1.95*nx, ty + sgn*1.95*ny
        star = star.union(bar(ex, ey, qx, qy, 2.4, 2.0, zs))
        peg = cq.Workplane("XY").circle(0.75).extrude(zd + 2.5 - (zs + 2.0)).translate((qx, qy, zs + 2.0))
        pegs = peg if pegs is None else pegs.union(peg)
    qx, qy = tlx - 1.95*thx, tly - 1.95*thy                            # peg against the detent L-tail
    star = star.union(bar(ex, ey, qx, qy, 2.4, 2.0, zs))
    pegs = pegs.union(cq.Workplane("XY").circle(0.75).extrude(zd + 2.5 - (zs + 2.0)).translate((qx, qy, zs + 2.0)))
    star = star.cut(cq.Workplane("XY").circle(1.05).extrude(4).translate((ex, ey, zs - 1)))
    parts["release star"] = star
    parts["release star pegs"] = pegs
    parts["release star axle"] = cq.Workplane("XY").circle(1.0).extrude(zs + 2.0 - BOTTOM_Z[1]).translate((ex, ey, BOTTOM_Z[1]))
    return parts

def balance_assembly():
    bx, by = BAL_XY; zc = LV["L7"]
    ro, ri, th = dp.BAL_R_OUT*1e3, dp.BAL_R_IN*1e3, dp.BAL_THICK*1e3
    ring = cq.Workplane("XY").circle(ro).circle(ri).extrude(th).translate((bx, by, zc-th/2))
    # 4 radial trim slots in the rim (M3 brass nuts slide 23..29 mm) -> regulation
    for i in range(4):
        a = 45 + 90*i
        slot = cq.Workplane("XY").box(dp.R_TRIM_MAX*1e3 - dp.R_TRIM_MIN*1e3 + 1.5, 3.3, th+1, centered=(False, True, True)).translate((dp.R_TRIM_MIN*1e3 - 0.75, 0, 0)).rotate((0,0,0),(0,0,1), a).translate((bx, by, zc))
        ring = ring.cut(slot)
    hub = cq.Workplane("XY").circle(dp.HUB_R*1e3).extrude(th).translate((bx, by, zc-th/2))
    spokes = None
    for i in range(4):
        a = 90*i
        sp = cq.Workplane("XY").box(ri - dp.HUB_R*1e3 + 1, 2.5, 2.5, centered=(False, True, True)).translate((dp.HUB_R*1e3 - 0.5, 0, 0)).rotate((0,0,0),(0,0,1), a).translate((bx, by, zc))
        spokes = sp if spokes is None else spokes.union(sp)
    balance = ring.union(hub).union(spokes)
    # band clamp slot through the hub (band in the XZ plane, along Z)
    balance = balance.cut(cq.Workplane("XY").box(dp.BAND_W*1e3+0.2, dp.BAND_T*1e3+0.2, th+1).translate((bx, by, zc)))
    trims = None
    for i in range(4):
        a = 45 + 90*i
        n = cq.Workplane("XY").polygon(6, 5.5).extrude(2.4).faces(">Z").workplane().hole(3.0).translate((dp.R_TRIM*1e3, 0, th/2)).rotate((0,0,0),(0,0,1), a).translate((bx, by, zc))
        trims = n if trims is None else trims.union(n)
    # impulse pin on the hub, toward the escape wheel
    ex, ey = S["S7"]; ang = math.atan2(ey-by, ex-bx)
    pin = cq.Workplane("XY").circle(0.6).extrude(th/2 + 9.5).translate((bx + (dp.HUB_R*1e3+1.2)*math.cos(ang), by + (dp.HUB_R*1e3+1.2)*math.sin(ang), zc-9.5))
    # taut band: hub to upper clamp (deck underside boss) and to lower clamp (tensioner on bottom plate)
    Lh = dp.BAND_LH*1e3
    band_up = cq.Workplane("XY").box(dp.BAND_W*1e3, dp.BAND_T*1e3, Lh + th/2 + 2, centered=(True, True, False)).translate((bx, by, zc))
    band_dn = cq.Workplane("XY").box(dp.BAND_W*1e3, dp.BAND_T*1e3, Lh + th/2 + 2, centered=(True, True, False)).translate((bx, by, zc - Lh - th/2 - 2))
    # rev 4: both band clamps adjustable (symmetric free-length trim L_h 18..26.5 mm, design_params.BAND_LH_MIN/MAX)
    z_uc = zc + th/2 + Lh - 4                       # upper clamp block bottom
    upper_clamp = cq.Workplane("XY").box(10, 8, 6, centered=(True, True, False)).translate((bx, by, z_uc))
    upper_clamp = upper_clamp.cut(cq.Workplane("XY").box(dp.BAND_W*1e3+0.2, dp.BAND_T*1e3+0.2, 7).translate((bx, by, z_uc + 3)))
    Lh_min, Lh_max = dp.BAND_LH_MIN*1e3, dp.BAND_LH_MAX*1e3
    # slotted hanger plate from the deck underside (3 mm steel), two M3 slots for the clamp screws
    upper_post = cq.Workplane("XY").box(12, 3, DECK_Z[1] - (zc + th/2 + Lh_min - 6), centered=(True, False, False)).translate((bx, by + 4, zc + th/2 + Lh_min - 6))
    upper_post = upper_post.union(cq.Workplane("XY").box(16, 8, 3, centered=(True, False, False)).translate((bx, by + 4, DECK_Z[1])))  # flange on the deck top
    for sx in (-3.0, 3.0):
        slot = cq.Workplane("XZ").slot2D(Lh_max - Lh_min + 3.2, 3.2, 90).extrude(-5).translate((bx + sx, by + 3.5, zc + th/2 + (Lh_min + Lh_max)/2 - 1))
        upper_post = upper_post.cut(slot)
    screws = None
    for sx in (-3.0, 3.0):
        s = cq.Workplane("XZ").circle(1.5).extrude(-9).translate((bx + sx, by - 1, z_uc + 3))          # M3 x 10 (along +y)
        s = s.union(cq.Workplane("XZ").circle(2.7).extrude(-2).translate((bx + sx, by + 7, z_uc + 3)))  # head
        screws = s if screws is None else screws.union(s)
    z_lc = zc - th/2 - Lh - 4                       # lower clamp block bottom (floats on the tensioner spring)
    lower_clamp = cq.Workplane("XY").box(10, 8, 4, centered=(True, True, False)).translate((bx, by, z_lc))
    lower_clamp = lower_clamp.cut(cq.Workplane("XY").box(dp.BAND_W*1e3+0.2, dp.BAND_T*1e3+0.2, 5).translate((bx, by, z_lc + 2)))
    # tensioner: steel coil spring between the lower clamp and a spring-seat block clamped on a slotted bracket
    SPRING_LEN = 6.5
    z_seat = z_lc - SPRING_LEN                      # seat block top
    seat = cq.Workplane("XY").box(10, 8, 3, centered=(True, True, False)).translate((bx, by, z_seat - 3))
    tens = None
    for i in range(4):
        ring = cq.Workplane("XY").circle(3.4).circle(2.6).extrude(0.8).translate((bx, by, z_seat + 0.2 + i*1.55))
        tens = ring if tens is None else tens.union(ring)
    z_b0 = BOTTOM_Z[0]; z_b1 = zc - th/2 - Lh_min - 4 + 8
    bracket = cq.Workplane("XY").box(12, 3, z_b1 - z_b0, centered=(True, False, False)).translate((bx, by + 4, z_b0))
    for sx in (-3.0, 3.0):   # slot covers seat screws (z_lc-8) and lower-clamp guide pins (z_lc+2) over L_h 18..26.5 mm
        slot = cq.Workplane("XZ").slot2D(22.0, 3.2, 90).extrude(-5).translate((bx + sx, by + 3.5, zc - th/2 - 4 - (Lh_min + Lh_max)/2 - 3.0))
        bracket = bracket.cut(slot)
    for sx in (-14.0, 14.0):   # feet on the bottom plate either side of the clearance pocket
        bracket = bracket.union(cq.Workplane("XY").box(4, 12, 3, centered=(True, True, False)).translate((bx + sx, by, BOTTOM_Z[1])))
    bracket = bracket.union(cq.Workplane("XY").box(32, 3, 3, centered=(True, False, False)).translate((bx, by + 4, BOTTOM_Z[1])))
    for sx in (-3.0, 3.0):
        s = cq.Workplane("XZ").circle(1.5).extrude(-9).translate((bx + sx, by - 1, z_seat - 1.5))     # seat screws
        s = s.union(cq.Workplane("XZ").circle(2.7).extrude(-2).translate((bx + sx, by + 7, z_seat - 1.5)))
        g = cq.Workplane("XZ").circle(1.5).extrude(-8).translate((bx + sx, by - 1, z_lc + 2))          # lower-clamp guide pins (slide in slots)
        screws = screws.union(s).union(g)
    guide = bracket
    return dict(balance=balance, trims=trims, impulse_pin=pin, band_upper=band_up, band_lower=band_dn,
                upper_clamp=upper_clamp, upper_post=upper_post, lower_clamp=lower_clamp, tensioner=tens,
                tensioner_seat=seat, tensioner_guide=guide, clamp_screws=screws)

# ------------------------------------------------------------------ assembly
def build(cassette="10 mL/h", fill=0.85, exploded=0.0, tripped=False, balance_angle=0.0, xc_override=None, stem_advance=0.0):
    asm = cq.Assembly(name="ChronoFlow")
    def add(obj, name, color, dz=0.0):
        if obj is None: return
        if exploded: obj = obj.translate((0, 0, dz*exploded))
        asm.add(obj, name=name, color=cq.Color(color))
    bp = bottom_plate().cut(cq.Workplane("XY").box(18, 16, 10).translate((BAL_XY[0], BAL_XY[1] + 2, BOTTOM_Z[0] + 3)))  # tensioner pocket
    add(bp, "bottom plate", COL["frame"], -40)
    deck = deck_plate().cut(cq.Workplane("XY").slot2D(2*FLAG_ROD_R + 3.0, 2*FLAG_ROD_R + 1.0, 0).extrude(10).translate((DROP_X - 0.5, RAIL_Y_TRIP, -8)))
    deck = deck.cut(cq.Workplane("XY").box(12.4, 3.4, 10, centered=(True, False, False)).translate((BAL_XY[0], BAL_XY[1] + 3.8, -8)))  # hanger-plate slot
    add(deck, "deck plate", COL["frame"], 0)
    add(pillars(), "pillars", COL["steel"], -20)
    add(feet(), "feet", COL["rubber"], -50)
    add(side_window_frame(), "balance window", "#a9d6e5", -20)
    add(rail_posts(), "rail posts", COL["print"], 30)
    add(rails(), "rails 6 mm", COL["steel"], 30)
    add(cradle(), "syringe cradle + flange latch", COL["print"], 40)
    add(barrel_clamp(), "barrel clamp", COL["print2"], 70)
    add(load_run_lever(), "LOAD/RUN interlock lever", COL["flag"], 70)
    bar, stopper, plunger, thumb_x = syringe(fill)
    xc = thumb_x + 3 + 10 + 2 + CAR_LEN/2      # pusher jaw face against the thumb press
    if xc_override is not None: xc = xc_override
    assert CAR_X_MIN - 1 <= xc <= CAR_X_MAX + 1, xc
    add(bar, "syringe barrel (BD 20 mL)", COL["syringe"], 60)
    add(stopper, "plunger stopper", COL["rubber"], 60)
    add(plunger, "plunger", COL["syringe"], 60)
    add(carriage(xc), "carriage + rack", COL["print"], 50)
    tm, t_spacer = trip_mechanism(xc, tripped, stem_advance=stem_advance)
    tcol = dict(pusher=COL["steel"], spacer=COL["print2"], trip_spring=COL["spring"], sear=COL["print2"], sear_pin=COL["steel"],
                sear_boss=COL["print"], sear_spring=COL["spring"], flag_rod=COL["flag"], flag_boss=COL["print"], flag_spring=COL["spring"],
                rocker=COL["print2"], rocker_boss=COL["print"], rocker_pin=COL["steel"], rocker_spring=COL["spring"],
                capture_pin=COL["steel"], capture_boss=COL["print"])
    tname = dict(pusher="pusher head + trip stem", spacer=f"trip preload spacer ({t_spacer:.1f} mm)", trip_spring="occlusion trip spring",
                 sear="trip sear", sear_pin="sear pivot pin", sear_boss="sear boss", sear_spring="sear return spring",
                 flag_rod="flag rod + occlusion flag", flag_boss="flag rod guide", flag_spring="flag spring",
                 rocker="capture rocker", rocker_boss="rocker boss", rocker_pin="rocker pivot pin", rocker_spring="rocker torsion spring",
                 capture_pin="stem capture pin", capture_boss="capture pin guide")
    for k, p in tm.items(): add(p, tname[k], tcol[k], 50)
    rl_name = dict(trip_rail="trip rail", rail_posts="trip rail posts", rail_cranks="trip rail cranks", rail_pins="trip rail pivot pins",
                   stop_linkage="drop rod + balance stop fork", reset_ramp="trip reset ramp")
    for k, p in trip_rail_linkage(xc, tripped).items():
        add(p, rl_name[k], COL["steel"] if k in ("stop_linkage", "rail_pins") else COL["print2"], 80)
    drum, coil, bracket, axle = spring_drum()
    add(drum, "spring drum", COL["print"], 20); add(coil, "constant-force spring coil", COL["spring"], 20)
    add(bracket, "drum bracket", COL["print"], 20); add(axle, "drum axle", COL["steel"], 20)
    add(spring_tape(xc), "spring tape", COL["spring"], 20)
    for k, s in shafts().items(): add(s, f"shaft {k}", COL["steel"], 0)
    for k, g in gear_train(cassette, xc).items():
        add(g, k, COL["print2"] if "cassette" in k else (COL["brass"] if "escape" in k else COL["print"]), -10 if "cassette" in k else 0)
    add(ratchet_click(), "winding ratchet click", COL["steel"], 0)
    add(cassette_plate(), f"rate cassette plate ({cassette})", COL["print2"], -10)
    for k, p in spring_arm_escapement().items():
        add(p, k, COL["steel"] if ("pin" in k or "spring" in k or "link" in k or "post" in k) else COL["print2"], -30)
    for k, p in balance_assembly().items():
        c = {"balance": COL["print"], "trims": COL["brass"], "impulse_pin": COL["steel"], "band_upper": COL["band"], "band_lower": COL["band"],
             "upper_clamp": COL["frame"], "upper_post": COL["frame"], "lower_clamp": COL["frame"], "tensioner": COL["spring"],
             "tensioner_seat": COL["frame"], "tensioner_guide": COL["frame"], "clamp_screws": COL["steel"]}[k]
        if balance_angle and k in ("balance", "trims", "impulse_pin"):
            p = p.rotate((BAL_XY[0], BAL_XY[1], 0), (BAL_XY[0], BAL_XY[1], 1), math.degrees(balance_angle))
        add(p, f"balance {k}", c, -30)
    return asm, xc

def _safe(n):
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", n).strip("_")

def export_stls(asm, out_dir, manifest, tol=0.05, ang=0.2):
    """One STL per assembly child + a JSON manifest (name, file, colour); returns the names of non-watertight meshes."""
    import json, os
    os.makedirs(out_dir, exist_ok=True)
    man, bad = [], []
    for ch in asm.children:
        shp = ch.obj.val() if hasattr(ch.obj, "val") else ch.obj
        fn = f"{out_dir}/{_safe(ch.name)}.stl"
        cq.exporters.export(cq.Workplane().add(shp), fn, tolerance=tol, angularTolerance=ang)
        col = ch.color.toTuple()[:3] if ch.color is not None else (0.6, 0.6, 0.6)
        man.append(dict(name=ch.name, file=fn, color=list(col)))
        try:
            import trimesh
            if not trimesh.load(fn, force="mesh").is_watertight: bad.append(ch.name)
        except ImportError:
            pass
    json.dump(man, open(manifest, "w"), indent=1)
    return man, bad

def export_all(cassette="10 mL/h", fill=0.85):
    """STEP + GLB + per-part STLs (armed pose), STLs of the tripped pose for the figure, layout.json, STL zip."""
    import json, zipfile, shutil
    asm, xc = build(cassette=cassette, fill=fill)
    asm.save(f"{OUT}/chronoflow_assembly.step")
    asm.save(f"{OUT}/chronoflow_assembly.glb", tolerance=0.3, angularTolerance=0.5)
    shutil.rmtree(f"{OUT}/stl", ignore_errors=True); shutil.rmtree(f"{OUT}/stl_tripped", ignore_errors=True)
    man, bad = export_stls(asm, f"{OUT}/stl", f"{OUT}/manifest.json")
    asm_t, _ = build(cassette=cassette, fill=fill, tripped=True)
    export_stls(asm_t, f"{OUT}/stl_tripped", f"{OUT}/manifest_tripped.json")
    json.dump(dict(shafts=S, levels=LV, balance_xy=BAL_XY, deck_z=DECK_Z, bottom_z=BOTTOM_Z, carriage_xc=xc, n_parts=len(asm.children),
                   trip=dict(sense_travel_mm=SENSE_TRAVEL, flag_travel_mm=FLAG_TRAVEL, sear_ratio=SEAR_RATIO, crank_posts=CRANK_POSTS,
                             stop_angles_deg=STOP_ANGLES, stop_r_mm=STOP_R)), open(f"{OUT}/layout.json", "w"), indent=1)
    with zipfile.ZipFile(f"{OUT}/chronoflow_stl_parts.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for m in man: z.write(m["file"], arcname=m["file"].split("/")[-1])
        z.write(f"{OUT}/manifest.json", arcname="manifest.json")
    return asm, xc, man, bad

if __name__ == "__main__":
    asm, xc, man, bad = export_all()
    print("carriage centre", round(xc, 2), "children", len(asm.children), "non-watertight:", bad)
