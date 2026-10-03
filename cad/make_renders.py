"""Render the CAD panels of the manuscript (Fig. 1a overview, Fig. 6a-e) from the exported STL manifests.
Run after cad_model.export_all():  python make_renders.py   (needs trimesh + matplotlib; no CadQuery needed)."""
import json, os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd())
import render_cad as rc

PX = json.load(open("panels/panel_px.json"))
def fs(f, L, dpi=300): w, h = PX[f][L]; return (w/dpi, h/dpi)
LAY = json.load(open("cad_out/layout.json")); xc_ = LAY["carriage_xc"]; Sx = LAY["shafts"]; L = LAY["levels"]; bx, by = LAY["balance_xy"]
ex = ("balance window",)
parts = rc.load_parts("cad_out/manifest.json")
parts_tripped = rc.load_parts("cad_out/manifest_tripped.json")

# ---- Fig 1a: labelled overview (armed)
labels1 = [("constant-force spring drum (44.5 N)", (18, 45, 21), (0.2, 0.27)),
           ("syringe cradle + barrel clamp", (60, 45, 70), (0.17, 0.92)),
           ("LOAD/RUN interlock", (103, 67, 85), (0.36, 0.95)),
           ("BD 20 mL syringe", (90, 45, 66), (0.46, 0.82)),
           ("snap-action occlusion trip:\nsear, flag rod, capture pin", (xc_ + 13, 64, 84), (0.72, 0.92)),
           ("trip rail on parallel cranks\n(lifts the balance stop fork)", (xc_ - 40, 72, 102), (0.88, 0.8)),
           ("carriage + rack on 6 mm rails", (xc_, 70, 50), (0.88, 0.64)),
           ("rack pinion + winding ratchet", (158.0, 12.0, 32), (0.88, 0.5)),
           ("step-up train 5×5×5×2.5", (158.0, Sx["S3"][1], 8), (0.88, 0.32)),
           ("movement below deck:\ncassette · escapement · balance", (82.0, -62.0, -36.0), (0.3, 0.08))]
rc.render(parts, "panels/fig1_a.png", elev=30, azim=35, exclude=ex, figsize=fs("fig1", "a"), dpi=300, labels=labels1, exact=True, label_fs=6.5, pad=30)

# ---- Fig 6a: exploded view
offsets = {"deck plate": 0, "bottom plate": -40, "pillars": -20, "feet": -50, "rail posts": 45, "rails": 45, "cradle": 60, "barrel clamp": 100,
           "LOAD/RUN": 100, "syringe barrel": 85, "stopper": 85, "plunger": 85, "carriage": 55, "pusher": 55, "trip preload": 55, "trip spring": 55,
           "occlusion trip spring": 55, "sear": 60, "flag": 60, "rocker": 60, "capture": 60, "trip rail": 115, "drop rod": 115,
           "spring drum": 30, "spring coil": 30, "drum bracket": 30, "drum axle": 30, "spring tape": 30,
           "W1": 20, "P1": 20, "W2": 14, "P2": 14, "W3": 8, "P3": 8, "W4": -8, "P4": -8, "cassette": -30, "escape wheel": -55,
           "impulse arm": -60, "latch": -60, "wheel detent": -60, "balance": -70, "winding ratchet": 24, "P0": 28}
parts_ex = rc.load_parts("cad_out/manifest.json", offsets=offsets)
rc.render(parts_ex, "panels/fig5_a.png", elev=25, azim=35, exclude=ex, figsize=fs("fig5", "a"), dpi=300, exact=True, pad=8)

# ---- Fig 6b: escapement module (lower movement) with the wheel detent
lower_keys = ("shaft S4", "shaft S5", "shaft S6", "shaft S7", "cassette", "escape wheel", "impulse arm", "balance ", "wheel detent", "drop rod")
lower = [p for p in parts if (any(k in p["name"] for k in lower_keys) or p["name"] == "latch") and "window" not in p["name"]]
labels = [("escape wheel (15 T)", (Sx["S7"][0], Sx["S7"][1], L["L7"]), (0.86, 0.6)),
          ("rate cassette (4 gears)", (Sx["S6"][0], Sx["S6"][1], L["L5"]), (0.55, 0.95)),
          ("spring-arm impulse levers", (Sx["S7"][0] - 14, Sx["S7"][1] - 8, -43.0), (0.78, 0.18)),
          ("compliant wheel detent\n(k ≈ 2000 N/m flexure)", (Sx["S7"][0] + 20, Sx["S7"][1] - 2, -38.0), (0.92, 0.84)),
          ("balance stop fork\n(occlusion trip)", (bx + 3, by + 17.7, -41.0), (0.5, 0.08)),
          ("taut-band balance + trim nuts", (bx, by + 28, L["L7"]), (0.22, 0.86)),
          ("steel ribbon 38.1 µm × 2 mm", (bx, by, -64.0), (0.18, 0.3)),
          ("slotted clamp hanger\n(L_h 18–26.5 mm)", (bx, by + 5.5, -14.0), (0.2, 0.6))]
rc.render(lower, "panels/fig5_b.png", elev=35, azim=20, figsize=fs("fig5", "b"), dpi=300, labels=labels, exact=True, label_fs=6.0, pad=25)

# ---- Fig 6c/d: occlusion trip close-up, armed vs tripped (carriage region only)
trip_keys = ("carriage", "pusher", "trip preload", "occlusion trip spring", "sear", "flag", "rocker", "capture", "trip rail", "rails 6 mm", "plunger")
def trip_view(pp, fname, letter):
    sel = [p for p in pp if any(k in p["name"] for k in trip_keys) and "posts" not in p["name"] and "cranks" not in p["name"] and "pivot pins" not in p["name"]]
    # crop the long rails/rail/plunger to the carriage neighbourhood
    import numpy as np
    out = []
    for p in sel:
        V, F = p["V"], p["F"]
        if any(k in p["name"] for k in ("rails 6 mm", "trip rail", "plunger")):
            cen = V[F].mean(axis=1); keep = (cen[:, 0] > xc_ - 45) & (cen[:, 0] < xc_ + 55)
            F = F[keep]
            if len(F) == 0: continue
        out.append(dict(name=p["name"], V=V, F=F, color=p["color"]))
    labs = [("sear (vertical pivot)", (xc_ + 20, 51, 85), (0.76, 0.93)),
            ("flag rod + 2 N spring", (xc_ + 13, 67, 90), (0.78, 0.76)),
            ("trip pin on stem collar", (xc_ + 44.6, 49, 80), (0.26, 0.93)),
            ("preload spacer (keyed)", (xc_ + 22, 45, 48), (0.22, 0.22)),
            ("rocker → capture pin", (xc_ + 6, 50, 76), (0.78, 0.38)),
            ("Ø8 stem + capture slot", (xc_ + 6, 45, 59), (0.5, 0.07))]
    rc.render(out, fname, elev=28, azim=125, figsize=fs("fig5", letter), dpi=300, labels=labs if letter == "c" else
              [("stem retracted 0.35 mm,\nsear turned, rod up 6 mm,\npin in slot → stem captured", (xc_ + 6, 45, 60), (0.5, 0.08))],
              exact=True, label_fs=6.0, pad=6)
trip_view(parts, "panels/fig5_c.png", "c")
trip_view(parts_tripped, "panels/fig5_d.png", "d")

# ---- Fig 6e: adjustable band clamps (balance column only)
clamp = [p for p in parts if (p["name"].startswith("balance ") and "window" not in p["name"]) or "bottom plate" in p["name"] or "deck plate" in p["name"]]
import numpy as np
cl = []
for p in clamp:
    V, F = p["V"], p["F"]
    if "plate" in p["name"]:
        cen = V[F].mean(axis=1); keep = (abs(cen[:, 0] - bx) < 18) & (abs(cen[:, 1] - by) < 18); F = F[keep]
        if len(F) == 0: continue
    cl.append(dict(name=p["name"], V=V, F=F, color=p["color"]))
labs = [("hanger plate with 2 slots\n(upper clamp, M3 × 2)", (bx + 6, by + 5.5, -12), (0.72, 0.9)),
        ("upper clamp block", (bx, by - 4, -12.5), (0.22, 0.78)),
        ("taut band", (bx, by, -40), (0.2, 0.55)),
        ("lower clamp (guide pins\nslide in the bracket slots)", (bx, by - 4, -62.5), (0.25, 0.3)),
        ("tensioner spring + seat block\n(M3 × 2 in slots)", (bx, by, -70), (0.72, 0.2)),
        ("pocket in bottom plate", (bx + 8, by + 6, -77), (0.72, 0.06))]
rc.render(cl, "panels/fig5_e.png", elev=12, azim=-55, figsize=fs("fig5", "e"), dpi=300, labels=labs, exact=True, label_fs=6.0, pad=6)
print("renders done")
