"""Orthographic painter's-algorithm renderer for the ChronoFlow part meshes (no OpenGL needed)."""
import json, math
import numpy as np, trimesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def load_parts(manifest="cad_out/manifest.json", offsets=None):
    parts = []
    for m in json.load(open(manifest)):
        mesh = trimesh.load(m["file"], force="mesh")
        if len(mesh.faces) < 5000:
            try:
                mesh = mesh.subdivide_to_size(max_edge=6.0, max_iter=6)
            except Exception:
                pass
        v = mesh.vertices.copy()
        if offsets:
            for key, dz in offsets.items():
                if key in m["name"]:
                    v = v + np.array([0, 0, dz])
                    break
        parts.append(dict(name=m["name"], V=v, F=mesh.faces, color=np.array(m["color"])))
    return parts

def rotation(elev, azim):
    e, a = math.radians(elev), math.radians(azim)
    Rz = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
    t = e - math.pi / 2
    Rx = np.array([[1, 0, 0], [0, math.cos(t), -math.sin(t)], [0, math.sin(t), math.cos(t)]])
    return Rx @ Rz   # world -> camera: x right, y up (world Z up on screen), z toward viewer

def render(parts, fname, elev=30, azim=-50, exclude=(), include=None, figsize=(9, 6.5), dpi=220,
           light=(0.3, 0.5, 0.8), title=None, labels=None, clip_y_max=None, clip_z_min=None, edge=False, exact=False, label_fs=8, pad=5):
    R = rotation(elev, azim)
    L = np.array(light); L = L / np.linalg.norm(L)
    polys, cols, depth = [], [], []
    for p in parts:
        if any(x in p["name"] for x in exclude): continue
        if include is not None and not any(x in p["name"] for x in include): continue
        V, F = p["V"], p["F"]
        if clip_y_max is not None or clip_z_min is not None:
            keep = np.ones(len(F), bool)
            cen = V[F].mean(axis=1)
            if clip_y_max is not None: keep &= cen[:, 1] <= clip_y_max
            if clip_z_min is not None: keep &= cen[:, 2] >= clip_z_min
            F = F[keep]
            if len(F) == 0: continue
        tri = V[F]                                  # (n,3,3)
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        nn = np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
        n = n / nn
        cam = tri @ R.T                              # (n,3,3)
        ncam = n @ R.T
        front = ncam[:, 2] > -0.2                     # cull strongly back-facing
        shade = 0.35 + 0.65 * np.clip(np.abs(n @ L), 0, 1)
        c = np.clip(p["color"][None, :] * shade[:, None], 0, 1)
        polys.append(cam[front][:, :, :2]); cols.append(c[front]); depth.append(cam[front][:, :, 2].mean(axis=1))
    P = np.concatenate(polys); C = np.concatenate(cols); D = np.concatenate(depth)
    order = np.argsort(D)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    pc = PolyCollection(P[order], facecolors=C[order], edgecolors=C[order], linewidths=0.25, antialiased=True)
    ax.add_collection(pc)
    ax.set_xlim(P[:, :, 0].min() - pad, P[:, :, 0].max() + pad); ax.set_ylim(P[:, :, 1].min() - pad, P[:, :, 1].max() + pad)
    ax.set_aspect("equal"); ax.axis("off")
    if title: ax.set_title(title, loc="left", fontsize=11)
    if labels:
        for text, (wx, wy, wz), (tx, ty) in labels:
            cx, cy, _ = R @ np.array([wx, wy, wz])
            tc = "axes fraction" if (abs(tx) <= 1.0 and abs(ty) <= 1.0) else "offset points"
            ax.annotate(text, xy=(cx, cy), xytext=(tx, ty), textcoords=tc, fontsize=label_fs,
                        arrowprops=dict(arrowstyle="-", lw=0.6, color="#222222"), ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888888", lw=0.5))
    if exact:
        fig.subplots_adjust(0, 0, 1, 1)
        fig.savefig(fname, facecolor="white")
    else:
        fig.savefig(fname, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return fname
