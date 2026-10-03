"""Pairwise solid-interference check of a CadQuery Assembly (bounding-box prefilter + OCC common)."""
import itertools, json
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

# intended contacts: press-fit pins/shafts in bores, band in clamps, springs seated on guides, trims on screws
INTENDED = [("shaft", None), ("pin", None), ("band", "clamp"), ("band", "balance"), ("tensioner", None), ("axle", None),
            ("rails", "posts"), ("rails", "carriage"), ("rails", "cradle"), ("trims", "balance"), ("stem", "carriage"),
            ("stopper", "barrel"), ("plunger", "stopper"), ("coil", "drum"), ("tape", "coil"), ("tape", "carriage"),
            ("pallet pin", "escape wheel"), ("arm", "escape wheel"), ("flexure", "arm"), ("click", "W1")]

def is_intended(a, b):
    for k1, k2 in INTENDED:
        if k2 is None:
            if k1 in a or k1 in b: return True
        else:
            if (k1 in a and k2 in b) or (k1 in b and k2 in a): return True
    return False

def check(asm, vol_tol=0.5, bb_tol=0.05):
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
            hits.append(dict(a=objs[i][0], b=objs[j][0], volume_mm3=round(v, 2), intended=is_intended(objs[i][0], objs[j][0])))
    return hits
