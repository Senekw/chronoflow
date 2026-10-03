"""Pairwise solid-interference check of the ChronoFlow CadQuery assembly (bounding-box prefilter + OCC common),
run over a set of mechanism poses.

rev 4 (review round 5): the earlier keyword whitelist ("shaft", "pin", "axle", "tensioner" matched ANY pair) is replaced
by an explicit list of allowed part pairs, each with the physical reason for the modelled overlap and a volume cap.
Everything else that overlaps by more than VOL_TOL mm^3 is reported as a fault.
"""
import itertools, json, math, sys
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

VOL_TOL = 0.05   # mm^3 below which a boolean common is treated as numerical noise / line contact

# (part A substring, part B substring, max allowed mm^3, reason).  Substrings must identify the pair uniquely.
ALLOWED = [
    ("impulse arm A", "impulse arm A pivot pin", 8.0, "pivot pin pressed into the arm root (flexure-mounted arm)"),
    ("impulse arm B", "impulse arm B pivot pin", 8.0, "pivot pin pressed into the arm root (flexure-mounted arm)"),
    ("impulse arm A", "impulse arm A pallet pin", 1.0, "pallet pin pressed into the arm"),
    ("impulse arm B", "impulse arm B pallet pin", 1.0, "pallet pin pressed into the arm"),
    ("escape wheel 15T", "impulse arm A pallet pin", 0.7, "locked tooth resting on the pallet pin (line contact drawn with 0.1 mm penetration)"),
    ("escape wheel 15T", "impulse arm B pallet pin", 0.7, "locked tooth resting on the pallet pin"),
    ("balance balance", "balance impulse_pin", 3.0, "impulse pin pressed into the hub"),
    ("wheel detent", "wheel detent pin", 1.0, "locking pin pressed into the detent lever"),
    ("wheel detent flexure", "wheel detent post", 7.0, "flexure root anchored in its post (printed as one part)"),
    ("balance upper_clamp", "balance clamp_screws", 75.0, "M3 screws threaded through the clamp block"),
    ("balance lower_clamp", "balance clamp_screws", 75.0, "guide pins pressed into the lower clamp"),
    ("balance tensioner_seat", "balance clamp_screws", 75.0, "M3 screws threaded through the seat block"),
    ("balance tensioner_guide", "balance clamp_screws", 45.0, "screw shanks passing through the bracket slots (slot drawn to the nominal position)"),
    ("balance balance", "balance band_upper", 0.2, "band clamped in the hub slot: the straight band is the untwisted reference; in the swung poses the band twists (not modelled)"),
    ("balance balance", "balance band_lower", 0.2, "as above"),
    ("carriage + rack", "P0 rack pinion", 1.0, "trapezoidal rack vs involute pinion tooth-tip overlap at the drawn phase (<= 0.11 mm^2 section)"),
    ("W1 60T", "P1 12T", 0.7, "tooth-tip overlap of the unshifted 12T pinion (profile shift +0.3 recommended)"),
    ("W2 60T", "P2 12T", 0.4, "tooth-tip overlap of the unshifted 12T pinion"),
    ("W3 60T", "P3 12T", 0.4, "tooth-tip overlap of the unshifted 12T pinion"),
    ("W4 30T", "P4 12T", 0.4, "tooth-tip overlap of the unshifted 12T pinion"),
]

def allowed_for(a, b):
    for ka, kb, cap, why in ALLOWED:
        if (ka in a and kb in b) or (ka in b and kb in a):
            return cap, why
    return None, None

def check(asm, vol_tol=VOL_TOL, bb_tol=0.05):
    objs = [(ch.name, ch.obj.val() if hasattr(ch.obj, "val") else ch.obj) for ch in asm.children]
    bbs = [o.BoundingBox() for _, o in objs]
    def ov(a, b):
        return not (a.xmax < b.xmin - bb_tol or b.xmax < a.xmin - bb_tol or a.ymax < b.ymin - bb_tol or
                    b.ymax < a.ymin - bb_tol or a.zmax < b.zmin - bb_tol or b.zmax < a.zmin - bb_tol)
    hits = []
    for i, j in itertools.combinations(range(len(objs)), 2):
        if not ov(bbs[i], bbs[j]): continue
        common = BRepAlgoAPI_Common(objs[i][1].wrapped, objs[j][1].wrapped); common.Build()
        if not common.IsDone(): continue
        props = GProp_GProps(); BRepGProp.VolumeProperties_s(common.Shape(), props)
        v = props.Mass()
        if v > vol_tol:
            cap, why = allowed_for(objs[i][0], objs[j][0])
            ok = cap is not None and v <= cap
            hits.append(dict(a=objs[i][0], b=objs[j][0], volume_mm3=round(v, 3), allowed=ok,
                             reason=why if cap is not None else "UNEXPECTED", cap_mm3=cap))
    return hits

def poses():
    """Carriage full/half/empty x armed/tripped x balance at rest and at +/- amplitude (0.6 rad)."""
    P = []
    for fill in (1.0, 0.5, 0.0):
        for ang in (0.0, 0.6, -0.6):
            P.append(dict(fill=fill, tripped=False, balance_angle=ang))
        P.append(dict(fill=fill, tripped=True, balance_angle=0.0))
    return P

def run(cassette="10 mL/h", out="results/interference_report.json"):
    import cad_model as cm
    report = []
    for p in poses():
        asm, xc = cm.build(cassette=cassette, **p)
        hits = check(asm)
        faults = [h for h in hits if not h["allowed"]]
        report.append(dict(pose=p, xc=round(xc, 2), n_parts=len(asm.children), n_contacts=len(hits), faults=faults,
                           allowed=[h for h in hits if h["allowed"]]))
        print(f"fill {p['fill']:.2f} xc {xc:6.2f} tripped {p['tripped']!s:5} bal {p['balance_angle']:+.1f}: "
              f"{len(hits)} contacts, {len(faults)} faults" + ("" if not faults else " -> " + "; ".join(f"{h['a']} | {h['b']} {h['volume_mm3']}" for h in faults)))
    json.dump(dict(cassette=cassette, vol_tol=VOL_TOL, allowed_pairs=[dict(a=a, b=b, cap_mm3=c, reason=r) for a, b, c, r in ALLOWED],
                   poses=report), open(out, "w"), indent=1)
    return report


def _common_volume(a, b):
    c = BRepAlgoAPI_Common(a.wrapped, b.wrapped); c.Build()
    props = GProp_GProps(); BRepGProp.VolumeProperties_s(c.Shape(), props)
    return props.Mass()

def verify_trip(fill=0.5, out="results/trip_capture_check.json"):
    """Functional checks of the modelled snap-action trip by solid booleans:
    (1) armed: capture pin clear of the stem; (2) tripped: pin engaged in the stem slot;
    (3) tripped + stem advanced 0.1/0.2 mm toward the syringe: pin/stem overlap (captured -> no release bolus);
    (4) tripped: balance free to swing +/-0.04 rad inside the fork but blocked at +/-0.08 rad (< 6 deg unlock angle)."""
    import cad_model as cm
    res = {}
    def parts(**kw):
        asm, xc = cm.build(fill=fill, **kw); return {c.name: c.obj.val() for c in asm.children}, xc
    P, xc = parts(tripped=False)
    res["armed_pin_stem_overlap_mm3"] = round(_common_volume(P["stem capture pin"], P["pusher head + trip stem"]), 4)
    pin_bb = P["stem capture pin"].BoundingBox(); res["armed_pin_tip_minus_stem_top_mm"] = round(pin_bb.zmin - (cm.SYR_Z + cm.STEM_R), 3)
    P, _ = parts(tripped=True)
    res["tripped_pin_stem_overlap_mm3"] = round(_common_volume(P["stem capture pin"], P["pusher head + trip stem"]), 4)
    pin_bb = P["stem capture pin"].BoundingBox(); res["tripped_pin_engagement_mm"] = round((cm.SYR_Z + cm.STEM_R) - pin_bb.zmin, 3)
    for adv in (0.1, 0.2):
        P, _ = parts(tripped=True, stem_advance=adv)
        res[f"tripped_stem_advanced_{adv}mm_overlap_mm3"] = round(_common_volume(P["stem capture pin"], P["pusher head + trip stem"]), 4)
    for ang in (0.04, 0.08, -0.04, -0.08):
        P, _ = parts(tripped=True, balance_angle=ang)
        res[f"tripped_balance_{ang:+.2f}rad_fork_overlap_mm3"] = round(_common_volume(P["drop rod + balance stop fork"], P["balance balance"]), 4)
    P, _ = parts(tripped=False, balance_angle=0.6)
    res["armed_balance_+0.60rad_fork_overlap_mm3"] = round(_common_volume(P["drop rod + balance stop fork"], P["balance balance"]), 4)
    res["sense_travel_mm"] = cm.SENSE_TRAVEL; res["flag_travel_mm"] = cm.FLAG_TRAVEL; res["xc"] = round(xc, 2)
    json.dump(res, open(out, "w"), indent=1)
    return res

if __name__ == "__main__":
    run(); print(verify_trip())
