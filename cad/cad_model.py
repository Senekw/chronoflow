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
    """Rack along X, teeth pointing -Y, pitch line at y=0, body from y=0 to +height."""
    p = math.pi*m; ha = m; hf = 1.25*m
    t = math.tan(pa)
    pts = [(-length/2, height)]
    x = -length/2
    while x < length/2 - 1e-6:
        # one pitch: root(-hf) ... tip(+ha) trapezoid, tooth centred at x + p/2
        x0 = x; xc = x + p/2
        pts += [(x0, -hf*0 + 0), ]  # placeholder replaced below
        pts.pop()
        pts += [(xc - p/4 - hf*t, hf), (xc - p/4 + ha*t, -ha), (xc + p/4 - ha*t, -ha), (xc + p/4 + hf*t, hf)]
        x += p
    pts += [(length/2, hf), (length/2, height)]
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
    lug = cq.Workplane("XY").box(10, 36, SYR_Z+18, centered=(False, True, False)).translate((98, SYR_Y, 0))
    for y in RAIL_Y:
        lug = lug.cut(cq.Workplane("YZ").circle(3.1).extrude(12).translate((97, y, RAIL_Z)))
    slot = cq.Workplane("XY").box(3.6, 25, 40, centered=(False, True, False)).translate((100, SYR_Y, SYR_Z-18))
    lug = lug.cut(slot)
    bore = cq.Workplane("YZ").circle(11.2).extrude(12).translate((97, SYR_Y, SYR_Z))
    lug = lug.cut(bore)
    lug = lug.cut(tape)
    return body.union(lug)

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
    rk = rack(110.0, face=RACK_Z[1]-RACK_Z[0], height=6.0)
    rk = rk.translate((xc, RACK_Y_PITCH, RACK_Z[0]))
    # rack body occupies y from pitch line to +8 -> attach to carriage underside (y 18..26)
    return body.union(rk)

def pusher_head(xc):
    """Pusher plate with plunger-flange gripper, mounted on the trip stem."""
    x_face = xc - CAR_LEN/2 - 2.0
    plate = cq.Workplane("YZ").rect(24, 32).extrude(4).translate((x_face-4, SYR_Y, SYR_Z))
    # gripper: a C-slot that captures the thumb-press disc (BD 20 mL ~ 3 mm thick)
    jaw = cq.Workplane("YZ").rect(24, 32).extrude(3).translate((x_face-10, SYR_Y, SYR_Z))
    jaw = jaw.cut(cq.Workplane("YZ").rect(18, 40).extrude(3).translate((x_face-10, SYR_Y, SYR_Z-10)))
    bridge = cq.Workplane("XY").box(6, 24, 4, centered=(False, True, False)).translate((x_face-10, SYR_Y, SYR_Z+12))
    stem = cq.Workplane("YZ").circle(4.0).extrude(CAR_LEN+14).translate((x_face-4, SYR_Y, SYR_Z))
    return plate.union(jaw).union(bridge).union(stem)

def trip_spring(xc):
    """Calibrated occlusion-trip compression spring behind the stem (helix as swept circle)."""
    x0 = xc + CAR_LEN/2 + 1
    # coil represented as a stack of 6 rings (render-friendly stand-in for the helix)
    sp = None
    for i in range(6):
        ring = cq.Workplane("YZ").circle(6.0).circle(4.8).extrude(0.9).translate((x0 + i*1.6, SYR_Y, SYR_Z))
        sp = ring if sp is None else sp.union(ring)
    return sp

def trip_knob_and_flag(xc):
    knob = cq.Workplane("YZ").polygon(6, 14).extrude(6).translate((xc + CAR_LEN/2 + 11, SYR_Y, SYR_Z))
    flag = cq.Workplane("XY").box(3, 8, 10, centered=(True, True, False)).translate((xc + CAR_LEN/2 - 6, SYR_Y, 66))
    return knob, flag

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
    fluid_len = dp.SYR_STROKE * fill_fraction          # full syringe: stopper near the flange end
    stopper_x = x_fl - 88 + 3 + fluid_len               # front face of the rubber stopper
    rod_len = 82.0                                      # empty syringe -> thumb press ~5 mm behind the flange
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

def gear_train(cassette):
    z1, z2, z3, z4 = dp.CASSETTES[cassette]
    G = {}
    f8 = dp.GEAR_FACE_STAGE1; f4 = dp.GEAR_FACE
    # rack pinion: rack tooth centred at x = S1.x (rack built centred on xc; phase so a pinion gap faces -Y... use 15 deg)
    G["P0 rack pinion"] = at(involute_gear(dp.Z_RACK_PINION, face=f8, phase_deg=(90 + 180/12) % 30), *S["S1"], RACK_Z[0])
    pairs = [("W1 60T", 60, "S1", "P1 12T", 12, "S2", f8, LV["L1"]), ("W2 60T", 60, "S2", "P2 12T", 12, "S3", f4, LV["L2"]),
             ("W3 60T", 60, "S3", "P3 12T", 12, "S4", f4, LV["L3"]), ("W4 30T", 30, "S4", "P4 12T", 12, "S5", f4, LV["L4"]),
             (f"cassette g1 {z1}T", z1, "S5", f"cassette g2 {z2}T", z2, "S6", f4, LV["L5"]),
             (f"cassette g3 {z3}T", z3, "S6", f"cassette g4 {z4}T", z4, "S7", f4, LV["L6"])]
    for nw, zw, sw, npn, zp, sp, f, lv in pairs:
        phw, php = mesh_phases(zw, zp, S[sw], S[sp])
        G[nw] = at(involute_gear(zw, face=f, phase_deg=phw), *S[sw], lv - f/2)
        G[npn] = at(involute_gear(zp, face=f, phase_deg=php), *S[sp], lv - f/2)
    G["escape wheel 15T"] = at(escape_wheel(face=f4), *S["S7"], LV["L7"]-f4/2)
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
    parts = {}
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
        # flexure spring: thin printed cantilever from a frame post pressing the arm toward the balance
        sx, sy = px + sgn*8*nx, py + sgn*8*ny
        parts[name + " flexure"] = cq.Workplane("XY").box(0.8, 14, 2.5, centered=(True, True, False)).rotate((0,0,0),(0,0,1), math.degrees(math.atan2(uy, ux))).translate((sx + 7*ux, sy + 7*uy, z0))
        parts[name + " post"] = cq.Workplane("XY").circle(2).extrude(16).translate((sx, sy, z0-8))
    # latch block between the arms, released by the balance impulse pin
    lx, ly = ex + 30*ux, ey + 30*uy
    parts["latch"] = cq.Workplane("XY").box(4, 10, 5.5, centered=(True, True, False)).rotate((0,0,0),(0,0,1), math.degrees(math.atan2(uy, ux))).translate((lx, ly, LV["L7"] - 8.5))
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
    upper_clamp = cq.Workplane("XY").box(10, 8, 6, centered=(True, True, False)).translate((bx, by, zc + th/2 + Lh - 4)) # boss on deck underside
    # upper clamp must reach the deck underside at z=-6: extend
    upper_post = cq.Workplane("XY").box(10, 8, DECK_Z[0] - (zc + th/2 + Lh), centered=(True, True, False)).translate((bx, by, zc + th/2 + Lh))
    lower_clamp = cq.Workplane("XY").box(10, 8, 4, centered=(True, True, False)).translate((bx, by, zc - th/2 - Lh - 4))
    # tensioner: steel coil spring between lower clamp and bottom plate
    zt0 = BOTTOM_Z[1]; zt1 = zc - th/2 - Lh - 4
    tens = None
    n = max(2, int((zt1 - zt0 - 3.2) / 1.6))
    for i in range(n):
        ring = cq.Workplane("XY").circle(3.4).circle(2.6).extrude(0.8).translate((bx, by, zt0 + 3.2 + i*1.6))
        tens = ring if tens is None else tens.union(ring)
    guide = cq.Workplane("XY").box(14, 12, 3, centered=(True, True, False)).translate((bx, by, zt0))
    return dict(balance=balance, trims=trims, impulse_pin=pin, band_upper=band_up, band_lower=band_dn,
                upper_clamp=upper_clamp.union(upper_post), lower_clamp=lower_clamp, tensioner=tens, tensioner_guide=guide)

# ------------------------------------------------------------------ assembly
def build(cassette="10 mL/h", fill=0.85, exploded=0.0):
    asm = cq.Assembly(name="ChronoFlow")
    def add(obj, name, color, dz=0.0):
        if obj is None: return
        if exploded: obj = obj.translate((0, 0, dz*exploded))
        asm.add(obj, name=name, color=cq.Color(color))
    add(deck_plate(), "deck plate", COL["frame"], 0)
    add(bottom_plate(), "bottom plate", COL["frame"], -40)
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
    assert CAR_X_MIN - 1 <= xc <= CAR_X_MAX + 1, xc
    add(bar, "syringe barrel (BD 20 mL)", COL["syringe"], 60)
    add(stopper, "plunger stopper", COL["rubber"], 60)
    add(plunger, "plunger", COL["syringe"], 60)
    add(carriage(xc), "carriage + rack", COL["print"], 50)
    add(pusher_head(xc), "pusher head + trip stem", COL["print2"], 50)
    add(trip_spring(xc), "occlusion trip spring", COL["spring"], 50)
    knob, flag = trip_knob_and_flag(xc)
    add(knob, "trip preload knob", COL["print2"], 50); add(flag, "occlusion flag", COL["flag"], 50)
    drum, coil, bracket, axle = spring_drum()
    add(drum, "spring drum", COL["print"], 20); add(coil, "constant-force spring coil", COL["spring"], 20)
    add(bracket, "drum bracket", COL["print"], 20); add(axle, "drum axle", COL["steel"], 20)
    add(spring_tape(xc), "spring tape", COL["spring"], 20)
    for k, s in shafts().items(): add(s, f"shaft {k}", COL["steel"], 0)
    for k, g in gear_train(cassette).items():
        add(g, k, COL["print2"] if "cassette" in k else (COL["brass"] if "escape" in k else COL["print"]), -10 if "cassette" in k else 0)
    add(ratchet_click(), "winding ratchet click", COL["steel"], 0)
    add(cassette_plate(), f"rate cassette plate ({cassette})", COL["print2"], -10)
    for k, p in spring_arm_escapement().items():
        add(p, k, COL["steel"] if "pin" in k else COL["print2"], -30)
    for k, p in balance_assembly().items():
        c = {"balance": COL["print"], "trims": COL["brass"], "impulse_pin": COL["steel"], "band_upper": COL["band"], "band_lower": COL["band"],
             "upper_clamp": COL["frame"], "lower_clamp": COL["frame"], "tensioner": COL["spring"], "tensioner_guide": COL["frame"]}[k]
        add(p, f"balance {k}", c, -30)
    return asm, xc

if __name__ == "__main__":
    asm, xc = build()
    asm.save(f"{OUT}/chronoflow_assembly.step")
    asm.save(f"{OUT}/chronoflow_assembly.glb", tolerance=0.3, angularTolerance=0.5)
    print("carriage centre", xc, "children", len(asm.children))
