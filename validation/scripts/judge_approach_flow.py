"""Approach flow of an empty-domain run of the wind-tunnel check: rows A1 to A3 of docs/VALIDATION.md, with the standard
error of every value from batch means of the averaging window, and the numbers the inflow calibration needs.

usage: uv run --with numpy python -I judge_approach_flow.py RUN_DIR GEOMETRY_DIR TUNNEL_JSON OUT_JSON [GAIN_USED|- [NEXT_GAIN|- [PREV_JUDGE_JSON|-]]] [--r R]
RUN_DIR/planes holds the x-normal planes of the run (xst-<frame>.f64, see report_approach_profiles.py; the method of the profiles is the
same: band mean over |y/b| <= 2, linear interpolation between layers to the 24 heights of the Inflow sheet, simulation divided
by the inflow speed at building height, tunnel by U_H). GAIN_USED is the calibration file the run's inflow table was made with.
NEXT_GAIN, when given, is written: the calibration file for the next run.

Rows (the pass marks are those of docs/VALIDATION.md, never changed here):
 A1  mean speed at the block position x/b = 0 against the sheet: within 5% from z/b = 0.125 up, within 10% at z/b = 0.0625.
 A2  turbulence intensity at building height at the block position, sigma_u / U at z/b = 2 of the simulation: 18% +/- 3 points.
 A3  the mean speed changes by less than 5% between x/b = -4.5 and x/b = -0.5 (the last two building heights of the fetch),
     at every height of the sheet.
Reported without a mark: the three sigmas at the block position against the sheet, their decay from the inlet plane x = 0
(the plane holds what the solver imposed), the Reynolds shear stress and its correlation coefficient at the block position.

Standard errors: every frame of the averaging window holds the sums so far; the difference of two
successive frames is one batch. A value is computed per batch exactly as for the whole window, and its standard error is the
standard deviation of the batch values over the root of their number. Batches must be much longer than the eddy time (7 s at
building height): 100 s batches are used. With fewer than three batches of equal length no error is given. The verdicts are
on the whole-window values; the errors are printed next to them.

Next calibration file:
 with --r R   mean-speed correction: gain "U" per sheet height = gain used x (tunnel U / simulated U
              at the block position) ** R, R between 0.5 and 1; the sigma gains are carried over unchanged.
 without      sigma correction: sigma gains = gains used x (sheet sigma / simulated sigma at the block position) ** (1 / S), S = 1
              or, with PREV_JUDGE_JSON, the measured sensitivity of the last two runs limited to 0.5 .. 1.5; "U" carried over.
 Gains are limited to 0.5 .. 2.5. Only the empty domain is read; no measuring point around the block enters.
"""
import json
import sys
from pathlib import Path

import numpy as np

argv = list(sys.argv[1:])
R_EXP = None
if "--r" in argv:
    i = argv.index("--r")
    R_EXP = float(argv[i + 1])
    del argv[i:i + 2]
    if not 0.5 <= R_EXP <= 1.0:
        raise SystemExit("--r must lie between 0.5 and 1")
RUN, GEOM, MEAS, OUT = (Path(a).resolve() for a in argv[:4])
opt = argv[4:] + ["-"] * 3
GAIN_USED = opt[0] if opt[0] != "-" else None
NEXT_GAIN = Path(opt[1]).resolve() if opt[1] != "-" else None
PREV = Path(opt[2]).resolve() if opt[2] != "-" else None
PL = RUN / "planes"
pm = dict(ln.split() for ln in (PL / "meta.txt").read_text().splitlines())
gm = json.loads((GEOM / "meta.json").read_text())
meas = json.loads(MEAS.read_text())
NB, cx, cy = gm["cells_per_b"], gm["centre_lattice"][0], gm["centre_lattice"][1]
u_ref, UH = float(pm["u_ref_lattice"]), meas["U_H_ms"]
sheet = sorted(meas["inflow"], key=lambda p: p["z_b"])
ZB = np.array([p["z_b"] for p in sheet])
T = {k: np.array([p[k] for p in sheet]) / UH for k in ("U", "su", "sv", "sw")}
Y_BAND = 2.0
j0, j1 = int(np.ceil(cy - Y_BAND * NB)), int(np.floor(cy + Y_BAND * NB)) + 1
KEYS = ("U", "su", "sv", "sw", "uw", "rho_uw")


def load(tag):
    meta = dict(ln.split(None, 1) for ln in (PL / ("xst_meta-%d.txt" % tag)).read_text().splitlines())
    planes = [int(v) for v in meta["planes"].split()]
    ny, nz, ns = int(meta["ny"]), int(meta["nz"]), int(meta["samples"])
    return {"tag": tag, "planes": planes, "n": ns, "S": np.memmap(PL / ("xst-%d.f64" % tag), np.float64, "r").reshape(len(planes), 7, nz, ny)}


def at_station(fr, x_b):
    """Sums at lattice x = cx + NB*x_b by linear interpolation between the two recorded planes; None when the frame lacks them."""
    xc = cx + NB * x_b
    lo = int(np.floor(xc + 1e-9))
    f = xc - lo
    need = [lo] if f < 1e-6 else [lo, lo + 1]
    if any(p not in fr["planes"] for p in need):
        return None
    a = np.array(fr["S"][fr["planes"].index(lo)])
    return a if f < 1e-6 else (1 - f) * a + f * np.array(fr["S"][fr["planes"].index(lo + 1)])


def at_heights(col):
    out = []
    for z_b in ZB:
        k = 0.5 + NB * z_b
        k0 = int(k)
        out.append(float(col[k0] * (1 - (k - k0)) + col[k0 + 1] * (k - k0)))
    return np.array(out)


def profile(S, n):
    """Band-mean profiles at the sheet heights, in U_H: U, the three sigmas, uw, and the correlation coefficient of u and w."""
    m = S[:3] / n
    var = np.maximum(S[3:6] / n - m ** 2, 0.0)
    uw = S[6] / n - m[0] * m[2]
    b = slice(j0, j1)
    U = at_heights(m[0][:, b].mean(axis=1)) / u_ref
    sig = [np.sqrt(at_heights(var[i][:, b].mean(axis=1))) / u_ref for i in range(3)]
    uwh = at_heights(uw[:, b].mean(axis=1)) / u_ref ** 2
    return {"U": U, "su": sig[0], "sv": sig[1], "sw": sig[2], "uw": uwh, "rho_uw": uwh / np.maximum(sig[0] * sig[2], 1e-12)}


def stderr(vals):
    """Standard error of the mean of the batch values (array over batches first); None with fewer than three batches."""
    v = np.asarray(vals, float)
    return None if v.shape[0] < 3 else v.std(axis=0, ddof=1) / np.sqrt(v.shape[0])


tags = sorted(int(p.stem.split("-")[1]) for p in PL.glob("xst_meta-*.txt"))
frames = sorted((load(t) for t in tags), key=lambda f: f["n"])
last = frames[-1]
ns = last["n"]
inlet = profile(at_station(last, -(cx / NB)), ns)                      # plane x = 0: what the solver imposed
STATIONS = (-10.0, -8.0, -6.0, -4.5, -4.0, -2.0, -0.5, 0.0)
st = {x: profile(at_station(last, x), ns) for x in STATIONS if at_station(last, x) is not None}
blk = st[0.0]

# batches: differences of successive frames, kept when their length equals the commonest one
cand = []
for prev, cur in zip([None] + frames[:-1], frames):
    n_b = cur["n"] - (prev["n"] if prev else 0)
    if n_b > 0:
        cand.append((prev, cur, n_b))
lens = [c[2] for c in cand]
n_batch = max(set(lens), key=lens.count) if lens else 0
batches = [c for c in cand if c[2] == n_batch]


def batch_profile(prev, cur, x_b):
    a = at_station(cur, x_b)
    if a is None:
        return None
    if prev is None:
        return profile(a, cur["n"])
    b = at_station(prev, x_b)
    return None if b is None else profile(a - b, cur["n"] - prev["n"])


bp = {x: [batch_profile(p, c, x) for p, c, _ in batches] for x in (0.0, -0.5, -4.5)}
ok = [k for k in range(len(batches)) if all(bp[x][k] is not None for x in bp)]
se = {k: stderr([bp[0.0][i][k] for i in ok]) for k in KEYS}
iH = int(np.argmin(np.abs(ZB - 2.0)))
ti = float(blk["su"][iH] / blk["U"][iH])
ti_se = stderr([bp[0.0][i]["su"][iH] / bp[0.0][i]["U"][iH] for i in ok])
chg = st[-0.5]["U"] / st[-4.5]["U"] - 1.0
chg_se = stderr([bp[-0.5][i]["U"] / bp[-4.5][i]["U"] - 1.0 for i in ok])


def f(v):
    return None if v is None else float(v)


def at(arr, i):
    return None if arr is None else float(arr[i])


dev = blk["U"] / T["U"] - 1.0
dev_se = None if se["U"] is None else se["U"] / T["U"]
up = ZB >= 0.125
iw = int(np.where(up)[0][np.abs(dev[up]).argmax()])
ilow = int(np.where(ZB == 0.0625)[0][0])
ic = int(np.abs(chg).argmax())
a1_up, a1_low = float(np.abs(dev[up]).max()), float(abs(dev[ilow]))
rows = {"A1": {"value": {"largest_deviation_from_z_b_0.125_up": a1_up, "at_z_b": float(ZB[iw]), "its_standard_error": at(dev_se, iw),
                         "deviation_at_z_b_0.0625": float(dev[ilow]), "standard_error_at_z_b_0.0625": at(dev_se, ilow)},
               "pass_mark": "within 5% from z/b = 0.125 up, within 10% at z/b = 0.0625", "pass": bool(a1_up <= 0.05 and a1_low <= 0.10)},
        "A2": {"value": ti, "standard_error": f(ti_se), "tunnel": float(T["su"][iH] / T["U"][iH]), "pass_mark": "18% +/- 3 points", "pass": bool(0.15 <= ti <= 0.21)},
        "A3": {"value": {"largest_change": float(np.abs(chg).max()), "at_z_b": float(ZB[ic]), "its_standard_error": at(chg_se, ic)},
               "pass_mark": "change below 5% between x/b = -4.5 and -0.5", "pass": bool(np.abs(chg).max() < 0.05)}}

gain_used = json.loads(Path(GAIN_USED).read_text()) if GAIN_USED else None
g_used = {k: (np.interp(ZB, gain_used["z_b"], gain_used[k]) if gain_used and k in gain_used else np.ones(len(ZB))) for k in ("su", "sv", "sw", "U")}
table = []
for i, z in enumerate(ZB):
    table.append({"z_b": float(z), "z_m": float(z * gm["b_m"]), "tunnel": {k: float(T[k][i]) for k in T}, "block": {k: float(blk[k][i]) for k in KEYS},
                  "block_standard_error": {k: at(se[k], i) for k in KEYS},
                  "inlet": {k: float(inlet[k][i]) for k in ("U", "su", "sv", "sw")}, "U_over_tunnel": float(blk["U"][i] / T["U"][i]), "U_over_tunnel_standard_error": at(dev_se, i),
                  "sigma_block_over_tunnel": {k: float(blk[k][i] / T[k][i]) for k in ("su", "sv", "sw")},
                  "sigma_block_over_inlet": {k: (float(blk[k][i] / inlet[k][i]) if inlet[k][i] > 1e-6 else None) for k in ("su", "sv", "sw")},
                  "sigma_inlet_over_target": {k: float(inlet[k][i] / (T[k][i] * g_used[k][i])) for k in ("su", "sv", "sw")},
                  "U_inlet_over_target": float(inlet["U"][i] / (T["U"][i] * g_used["U"][i])),
                  "U_change_last_two_heights": float(chg[i]), "U_change_standard_error": at(chg_se, i), "gain_used": {k: float(g_used[k][i]) for k in g_used}})
res = {"run": RUN.parent.name if RUN.name == "out" else RUN.name, "cells_per_b": NB, "samples": ns, "avg_start_s": float(pm["avg_start_s"]), "avg_end_s": float(pm["avg_end_s"]),
       "y_band_b": [-Y_BAND, Y_BAND], "gain_file_used": GAIN_USED,
       "batches": {"count": len(ok), "samples_each": n_batch, "frames": [c["tag"] for _, c, _ in batches], "note": "standard error = standard deviation of the batch values / root of their number"},
       "rows": rows, "table": table,
       "along_fetch": {("%g" % x): {k: st[x][k].tolist() for k in ("U", "su", "sv", "sw", "uw")} for x in st}, "z_b": ZB.tolist(),
       "playback_check": {"largest_deviation_of_inlet_sigma_from_target": {k: float(np.abs(inlet[k] / (T[k] * g_used[k]) - 1.0).max()) for k in ("su", "sv", "sw")},
                          "largest_deviation_of_inlet_U_from_target": float(np.abs(inlet["U"] / (T["U"] * g_used["U"]) - 1.0).max())}}

if NEXT_GAIN is not None:
    nxt = {"z_b": ZB.tolist()}
    if R_EXP is not None:      # mean-speed correction: the mean speed imposed at the inlet
        for k in ("su", "sv", "sw"):
            nxt[k] = g_used[k].tolist()
        nxt["U"] = np.clip(g_used["U"] * (T["U"] / blk["U"]) ** R_EXP, 0.5, 2.5).tolist()
        nxt["note"] = ("gain file for the run after %s: mean-speed gain U = gain used x (tunnel U / simulated U at the block position) ** %g per sheet height; "
                       "sigma gains unchanged; z_b = height over building width b" % (res["run"], R_EXP))
        nxt["r"] = R_EXP
    else:                      # sigma correction: the three sigmas
        sens = {k: np.ones(len(ZB)) for k in ("su", "sv", "sw")}
        if PREV is not None:
            prev = json.loads(PREV.read_text())
            for k in sens:
                g0 = np.array([r["gain_used"][k] for r in prev["table"]])
                s0 = np.array([r["block"][k] for r in prev["table"]])
                dg = np.log(g_used[k] / g0)
                moved = np.abs(dg) > 0.02
                sens[k][moved] = np.clip(np.log(blk[k][moved] / s0[moved]) / dg[moved], 0.5, 1.5)
        for k in ("su", "sv", "sw"):
            nxt[k] = np.clip(g_used[k] * (T[k] / blk[k]) ** (1.0 / sens[k]), 0.5, 2.5).tolist()
        if gain_used and "U" in gain_used:
            nxt["U"] = g_used["U"].tolist()
        nxt["note"] = "gains for the run after %s: gain used x (sheet / simulated at the block position) ** (1 / sensitivity); z_b = height over building width b" % res["run"]
        nxt["sensitivity"] = {k: sens[k].tolist() for k in sens}
    NEXT_GAIN.write_text(json.dumps(nxt, indent=1))
    res["next_gain_file"] = str(NEXT_GAIN.name)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(res, indent=1))


def pm_(v, e, fmt="%.3f", efmt="%.3f"):
    return (fmt % v) + (" ±" + (efmt % e) if e is not None else "      ")


print("%s: %d samples (%.0f-%.0f s), block position x/b = 0, band |y/b| <= %g; values in U_H; ± = standard error from %d batches of %d samples" %
      (res["run"], ns, res["avg_start_s"], res["avg_end_s"], Y_BAND, len(ok), n_batch))
print("  z/b    z[m] |  U: tunnel  sim ±se      dev ±se     | su: tunnel  sim ±se     ratio | sv: tunnel  sim ±se     ratio | sw: tunnel  sim ±se     ratio | block/inlet su sv sw | -uw     rho  | dU last 2H ±se")
for r in table:
    t, b, e = r["tunnel"], r["block"], r["block_standard_error"]
    des, ces = r["U_over_tunnel_standard_error"], r["U_change_standard_error"]
    print("  %6.4f %5.0f |   %.3f  %s %+5.1f%%%s |    %.3f  %s  %.2f |    %.3f  %s  %.2f |    %.3f  %s  %.2f |    %.2f %.2f %.2f    | %.5f %+.2f | %+5.1f%%%s" %
          (r["z_b"], r["z_m"], t["U"], pm_(b["U"], e["U"]), 100 * (r["U_over_tunnel"] - 1), (" ±%.1f" % (100 * des)) if des is not None else "     ",
           t["su"], pm_(b["su"], e["su"]), r["sigma_block_over_tunnel"]["su"], t["sv"], pm_(b["sv"], e["sv"]), r["sigma_block_over_tunnel"]["sv"],
           t["sw"], pm_(b["sw"], e["sw"]), r["sigma_block_over_tunnel"]["sw"],
           *(float("nan") if r["sigma_block_over_inlet"][k] is None else r["sigma_block_over_inlet"][k] for k in ("su", "sv", "sw")), -b["uw"], b["rho_uw"],
           100 * r["U_change_last_two_heights"], (" ±%.1f" % (100 * ces)) if ces is not None else ""))
pc = res["playback_check"]
print("inlet plane against its target (sheet x gain): largest deviation sigma_u %.1f%%, sigma_v %.1f%%, sigma_w %.1f%%; mean speed %.1f%%" %
      (100 * pc["largest_deviation_of_inlet_sigma_from_target"]["su"], 100 * pc["largest_deviation_of_inlet_sigma_from_target"]["sv"],
       100 * pc["largest_deviation_of_inlet_sigma_from_target"]["sw"], 100 * pc["largest_deviation_of_inlet_U_from_target"]))


def se_txt(v):
    return "no standard error (fewer than three equal batches)" if v is None else "standard error %.1f points" % (100 * v)


a = rows["A1"]["value"]
print("row A1: largest deviation %.1f%% (at z/b %g, %s) from z/b 0.125 up; %+.1f%% at z/b 0.0625 (%s)  -> %s" %
      (100 * a["largest_deviation_from_z_b_0.125_up"], a["at_z_b"], se_txt(a["its_standard_error"]), 100 * a["deviation_at_z_b_0.0625"], se_txt(a["standard_error_at_z_b_0.0625"]), "PASS" if rows["A1"]["pass"] else "FAIL"))
print("row A2: intensity at building height %.1f%% (%s; tunnel %.1f%%)  -> %s" % (100 * ti, se_txt(rows["A2"]["standard_error"]), 100 * rows["A2"]["tunnel"], "PASS" if rows["A2"]["pass"] else "FAIL"))
print("row A3: largest change of U between x/b -4.5 and -0.5: %.1f%% (at z/b %g, %s)  -> %s" %
      (100 * rows["A3"]["value"]["largest_change"], rows["A3"]["value"]["at_z_b"], se_txt(rows["A3"]["value"]["its_standard_error"]), "PASS" if rows["A3"]["pass"] else "FAIL"))
if NEXT_GAIN is not None:
    print("next gain file ->", NEXT_GAIN)
print("wrote", OUT)
