"""Approach-flow profiles from the x-normal planes of a run of the wind-tunnel check, against the tunnel's inflow sheet.

usage: uv run --with numpy python -I report_approach_profiles.py RUN_DIR GEOMETRY_DIR TUNNEL_JSON OUT_JSON [X_B ...]
RUN_DIR/planes holds meta.txt and, per output frame, xst-<frame>.f64 with xst_meta-<frame>.txt: for every recorded plane the
sums of u, v, w, uu, vv, ww, uw over the averaging window (lattice units, float64, ny*nz with y fastest). X_B are the
stations to report in building widths from the block centre (default -4 -2 -0.6 0); every other recorded station pair
is reported too.

Method. A station x/b lies at lattice x = cx + NB*x/b, between two recorded planes: the sums of the two planes are
interpolated linearly. Per cell: mean = sum/n, variance = sum of squares/n - mean^2, uw = sum(uw)/n - mean(u)*mean(w).
Then the mean over the band |y/b| <= 2 (the lateral extent of the measuring points; the side faces are 5.25 H away),
then linear interpolation in z to the 24 heights of the sheet (lattice z = 0.5 + NB*z/b, ground at z = 0.5).
sigma = sqrt(band mean of the variance). Normalisation: simulation by the inflow speed at building height (UREF),
tunnel by U(z = 2b) = U_H of the sheet. "imposed" is the same reading of the inlet plane x = 0, which holds the profile
the solver imposes. Settling: the same profiles from the last two quarters of the averaging window (differences of the
cumulative sums of the last three frames).
"""
import json
import sys
from pathlib import Path

import numpy as np

RUN, GEOM, MEAS, OUT = (Path(a).resolve() for a in sys.argv[1:5])
WANT = [float(a) for a in sys.argv[5:]] or [-4.0, -2.0, -0.6, 0.0]
PL = RUN / "planes"
pm = dict(ln.split() for ln in (PL / "meta.txt").read_text().splitlines())
gm = json.loads((GEOM / "meta.json").read_text())
meas = json.loads(MEAS.read_text())
NB, cx, cy = gm["cells_per_b"], gm["centre_lattice"][0], gm["centre_lattice"][1]
u_ref, UH = float(pm["u_ref_lattice"]), meas["U_H_ms"]
sheet = sorted(meas["inflow"], key=lambda p: p["z_b"])
ZB = [p["z_b"] for p in sheet]
Y_BAND = 2.0
j0, j1 = int(np.ceil(cy - Y_BAND * NB)), int(np.floor(cy + Y_BAND * NB)) + 1  # cells whose centres lie within |y/b| <= 2


def load(tag):
    meta = dict(ln.split(None, 1) for ln in (PL / ("xst_meta-%d.txt" % tag)).read_text().splitlines())
    planes = [int(v) for v in meta["planes"].split()]
    ny, nz, ns = int(meta["ny"]), int(meta["nz"]), int(meta["samples"])
    return planes, ns, np.fromfile(PL / ("xst-%d.f64" % tag), np.float64).reshape(len(planes), 7, nz, ny)


def at_station(planes, S, x_b):
    """Sums at lattice x = cx + NB*x_b by linear interpolation between the two recorded planes; None when not recorded."""
    xc = cx + NB * x_b
    lo = int(np.floor(xc + 1e-9))
    f = xc - lo
    if f < 1e-6:
        return (S[planes.index(lo)], [lo], 0.0) if lo in planes else None
    if lo in planes and lo + 1 in planes:
        return (1 - f) * S[planes.index(lo)] + f * S[planes.index(lo + 1)], [lo, lo + 1], f
    return None


def profile(S, n):
    """Band-mean statistics per layer z from the sums of one station: dict of arrays over z (lattice units)."""
    m = S[:3] / n
    var = np.maximum(S[3:6] / n - m ** 2, 0.0)
    uw = S[6] / n - m[0] * m[2]
    b = slice(j0, j1)
    return {"U": m[0][:, b].mean(axis=1), "V": m[1][:, b].mean(axis=1), "W": m[2][:, b].mean(axis=1), "var": var[:, :, b].mean(axis=2), "uw": uw[:, b].mean(axis=1),
            "U_spread_y": m[0][:, b].max(axis=1) - m[0][:, b].min(axis=1)}


def at_heights(col):
    """Linear interpolation of a per-layer column to the heights of the sheet."""
    out = []
    for z_b in ZB:
        k = 0.5 + NB * z_b
        k0 = int(k)
        out.append(float(col[k0] * (1 - (k - k0)) + col[k0 + 1] * (k - k0)))
    return np.array(out)


tags = sorted(int(p.stem.split("-")[1]) for p in PL.glob("xst_meta-*.txt"))
planes, ns, S = load(tags[-1])
imposed = at_heights(profile(S[planes.index(0)], ns)["U"]) / u_ref if 0 in planes else None
imp_col = profile(S[planes.index(0)], ns)["U"] / u_ref if 0 in planes else None
# every station the planes allow: the wanted ones first, then each recorded pair at its midpoint
cands = list(WANT) + [((a + 0.5) - cx) / NB for a in planes if a + 1 in planes and a > 2]
seen, stations = set(), []
for x_b in cands:
    key = round(x_b, 4)
    got = at_station(planes, S, x_b)
    if got is None or key in seen:
        continue
    seen.add(key)
    Sx, used, f = got
    p = profile(Sx, ns)
    U, su, sv, sw = at_heights(p["U"]) / u_ref, np.sqrt(at_heights(p["var"][0])) / u_ref, np.sqrt(at_heights(p["var"][1])) / u_ref, np.sqrt(at_heights(p["var"][2])) / u_ref
    uw = at_heights(p["uw"]) / u_ref ** 2
    rows = []
    for i, q in enumerate(sheet):
        rows.append({"z_b": q["z_b"], "z_m": q["z_b"] * gm["b_m"], "U": float(U[i]), "su": float(su[i]), "sv": float(sv[i]), "sw": float(sw[i]), "uw": float(uw[i]),
                     "U_over_imposed": float(U[i] / imposed[i]) if imposed is not None else None, "U_over_tunnel": float(U[i] / (q["U"] / UH))})
    stations.append({"x_b": x_b, "asked": x_b in WANT, "planes": used, "weight_of_second_plane": f, "profile": rows,
                     "first_cells": {"z_m": [(k - 0.5) * float(pm["dx_m"]) for k in (1, 2, 3, 4)], "U": [float(p["U"][k] / u_ref) for k in (1, 2, 3, 4)],
                                     "U_over_imposed": [float(p["U"][k] / u_ref / imp_col[k]) for k in (1, 2, 3, 4)] if imp_col is not None else None},
                     "U_spread_over_band_max": float((p["U_spread_y"][1:int(0.5 + NB * 7.5) + 2] / u_ref).max()),
                     "sigma_max": {"su": float(su.max()), "sv": float(sv.max()), "sw": float(sw.max())}, "minus_uw_max": float((-uw).max())})
stations.sort(key=lambda s: s["x_b"])

res = {"run": RUN.name, "dx_m": float(pm["dx_m"]), "cells_per_b": NB, "u_ref_lattice": u_ref, "samples": ns, "frame": tags[-1],
       "avg_start_s": float(pm["avg_start_s"]), "avg_end_s": float(pm["avg_end_s"]), "y_band_b": [-Y_BAND, Y_BAND], "y_cells": [j0, j1 - 1],
       "normalisation": "simulation by the inflow speed at building height (UREF), tunnel by U_H = %.3f m/s" % UH,
       "tunnel": [{"z_b": q["z_b"], "z_m": q["z_b"] * gm["b_m"], "U": q["U"] / UH, "su": q["su"] / UH, "sv": q["sv"] / UH, "sw": q["sw"] / UH} for q in sheet],
       "imposed_at_inlet": [{"z_b": z, "U": float(v)} for z, v in zip(ZB, imposed)] if imposed is not None else None, "stations": stations}

if len(tags) >= 3:  # settling between the last two quarters of the averaging window
    (_, n2, S2), (_, n1, S1) = load(tags[-2]), load(tags[-3])
    if ns > n2 > n1:
        d = []
        for st in stations:
            a = at_station(planes, S - S2, st["x_b"])[0]
            b = at_station(planes, S2 - S1, st["x_b"])[0]
            d.append(np.abs(at_heights(profile(a, ns - n2)["U"]) - at_heights(profile(b, n2 - n1)["U"])) / u_ref)
        d = np.array(d)
        res["settling"] = {"samples_last_quarter": ns - n2, "samples_quarter_before": n2 - n1, "mean_abs_dU": float(d.mean()), "max_abs_dU": float(d.max())}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(res, indent=1))
print("%s: %d samples (%.0f-%.0f s), band |y/b| <= %g, stations x/b: %s" % (RUN.name, ns, res["avg_start_s"], res["avg_end_s"], Y_BAND, ", ".join("%g" % s["x_b"] for s in stations)))
show = [s for s in stations if s["asked"]]
print("U/U_H   z/b    z(m)  tunnel imposed | " + "  ".join("x/b=%-5g" % s["x_b"] for s in show) + " || sigma_u/U_H: tunnel | " + "  ".join("x/b=%-5g" % s["x_b"] for s in show))
for i, q in enumerate(res["tunnel"]):
    print("      %6.4f %6.1f  %.3f  %s  | %s || %.3f | %s" % (q["z_b"], q["z_m"], q["U"], "%.3f" % imposed[i] if imposed is not None else "  -  ",
                                                                "  ".join("%.3f    " % s["profile"][i]["U"] for s in show), q["su"], "  ".join("%.4f   " % s["profile"][i]["su"] for s in show)))
for s in stations:
    print("x/b %6.2f: first cells (1, 3, 5, 7 m) U/U_H %s, over imposed %s; sigma max u %.4f v %.4f w %.4f; -uw max %.2e; spread of U across the band %.4f" %
          (s["x_b"], " ".join("%.3f" % v for v in s["first_cells"]["U"]), " ".join("%.3f" % v for v in (s["first_cells"]["U_over_imposed"] or [])),
           s["sigma_max"]["su"], s["sigma_max"]["sv"], s["sigma_max"]["sw"], s["minus_uw_max"], s["U_spread_over_band_max"]))
if "settling" in res:
    print("settling, last quarter against the one before: mean |dU| %.4f, max %.4f U_H" % (res["settling"]["mean_abs_dU"], res["settling"]["max_abs_dU"]))
print("wrote", OUT)
