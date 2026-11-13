"""Score a block run of the wind-tunnel check (AIJ case A): rows B1 to B6 of docs/VALIDATION.md on the 181 unique locations.

usage: uv run --with numpy python -I score_block_run.py RUN_DIR GEOMETRY_DIR TUNNEL_JSON OUT_JSON
RUN_DIR/planes holds meta.txt, vol_meta.txt and vol_mean.f32 (mean velocity in the 3D box over the averaging window, lattice
units, component-major, x fastest) and, when the run kept them, the means so far of earlier frames (vol_mean-<frame>.f32 with
vol_meta-<frame>.txt) and vol_std.f32 (standard deviations of u, v, w and the covariance of u and w; solver switch VOL_M2=1).
GEOMETRY_DIR holds meta.json and height.u16 of make_block_geometry.py; TUNNEL_JSON is the output of read_tunnel_data.py.

What is computed, and by which rule (the pass marks are those of docs/VALIDATION.md and are only read here):
 - Locations. The workbook has 186 points at 181 locations: five locations (y/b = 0, z/b = 0.125) were measured in two runs
   of the tunnel. Each location is scored once; its tunnel value is the mean of the two measurements (rule fixed before
   the first scored run; on the earlier version's result the scores are the same whichever of the two is kept). The pairs
   are reported as the tunnel's own repeatability.
 - Sampling: trilinear interpolation between the eight cell centres around the exact (x/b, y/b, z/b), as in the scorer of
   the earlier version (lattice x = cx + NB x/b, y = cy + NB y/b, z = 0.5 + NB z/b; solid cells count with zero velocity).
 - Normalisation: simulation by the inflow speed at building height (UREF), tunnel by U_H = U(z = 2b) of the Inflow sheet.
 - Hit (COST 732): |simulation - tunnel| <= 25% of |tunnel| or <= 0.05. Horizontal speed Uh = hypot(mean U, mean V); signed
   streamwise component U.
 - Regions: upstream x/b < -0.5; sides -0.5 <= x/b <= 0.5 and z/b < 2; roof -0.5 <= x/b <= 0.5 and
   z/b >= 2; wake below roof height x/b > 0.5 and z/b < 2; shear layer above the wake x/b > 0.5 and z/b >= 2; pedestrian
   level z/b = 0.125 (it overlaps the others).
 - Interval: Wilson score interval of a share of n points at 95% (z = 1.96). It treats the points as independent, which
   neighbouring points are not, so it is a lower bound of the real uncertainty.
 - FAC2: share of points with 0.5 <= simulation/tunnel <= 2 for the horizontal speed.
 - Reattachment: where the mean U on the centre line y/b = 0 behind the leeward face turns positive again, as in the diagnostics
   of the earlier version (simulation: the last negative cell column, linear interpolation to zero; tunnel: linear
   interpolation between the points of the vertical section), at z/b = 0.125, 0.5, 1.0, 1.5.
 - Settling (row B6): the horizontal speed at the 181 locations from two parts of the averaging window, mean and largest
   absolute difference. "quarters": the third and the fourth quarter of the averaging window (needs the means so far at 1/2
   and 3/4 of the samples); "halves": its two halves, which is what the earlier version's check reported as 0.015 and 0.091 (its
   window was half as long, so its halves were as long as the quarters here). Row B6 is judged on the quarters when the run
   kept them, else on the halves.
 - Fluctuation level (no pass mark): sigma_u, sigma_v, sigma_w at the locations (variance interpolated, then the root)
   against the tunnel's, as means per region and behind the block.
 - Standard errors (added 2026-11-10): the means so far of successive frames give batches of the
   averaging window (the difference of two frames). Per location: standard deviation of the batch values of U and of the
   horizontal speed over the root of the number of batches. Per hit rate and per reattachment length: delete-one-batch
   jackknife (the rate is recomputed from the mean of all batches but one; the batches are common to all locations, so the
   correlation between neighbouring points is kept). With fewer than three batches of equal length no error is given.
 - Row B6 is computed twice: as written (the two quarters), with its verdict; and by the rule that replaced it after the
   results were known (docs/VALIDATION.md), the standard error of the mean of the horizontal speed averaged over the
   locations against 0.025, printed without a verdict.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

Z95 = 1.959964
MARKS = {"B1": 0.66, "B2": 0.66, "B3": 0.60, "B4": 0.50, "B5": (0.99, 1.65), "B6": (0.01, 0.05)}  # pass marks of docs/VALIDATION.md; never edited here
REGIONS = ("upstream", "sides", "roof", "wake_below_roof", "shear_layer", "pedestrian")
B4_REGIONS = REGIONS


def region_of(x_b, z_b):
    if x_b < -0.5:
        return "upstream"
    if x_b > 0.5:
        return "wake_below_roof" if z_b < 2.0 else "shear_layer"
    return "sides" if z_b < 2.0 else "roof"


def wilson(hits, n):
    if n == 0:
        return None
    p, z2 = hits / n, Z95 * Z95
    centre, half = (p + z2 / (2 * n)) / (1 + z2 / n), Z95 * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return [max(0.0, centre - half), min(1.0, centre + half)]


def is_hit(pred, obs):
    return abs(pred - obs) <= 0.25 * abs(obs) or abs(pred - obs) <= 0.05


def unique_locations(points, u_h):
    """One record per location; the tunnel value of a location measured twice is the mean of the two measurements."""
    by_loc = {}
    for p in points:
        by_loc.setdefault((p["x_b"], p["y_b"], p["z_b"]), []).append(p)
    locs, pairs = [], []
    for (x_b, y_b, z_b), ps in by_loc.items():
        o = {k: float(np.mean([p[k] for p in ps])) / u_h for k in ("U", "V", "W", "su", "sv", "sw")}
        o["Uh"] = math.hypot(o["U"], o["V"])
        locs.append({"x_b": x_b, "y_b": y_b, "z_b": z_b, "sections": [p["section"] for p in ps], "region": region_of(x_b, z_b), "pedestrian": z_b == 0.125, "obs": o})
        if len(ps) > 1:
            pairs.append({"x_b": x_b, "y_b": y_b, "z_b": z_b, "sections": [p["section"] for p in ps],
                          "difference": {k: abs(ps[0][k] - ps[1][k]) / u_h for k in ("U", "V", "W", "su", "sv", "sw")}})
    return locs, pairs


class Box:
    """A field of the 3D box, normalised, with trilinear sampling at (x/b, y/b, z/b)."""

    def __init__(self, planes, name, meta_name, comps, scale, geom, hmap):
        vm = dict(ln.split(None, 1) for ln in (planes / meta_name).read_text().splitlines())
        self.x0, self.y0, self.nx, self.ny, self.nz, self.samples = (int(vm[k]) for k in ("x0", "y0", "nx", "ny", "nz", "samples"))
        self.a = np.fromfile(planes / name, np.float32).reshape(comps, self.nz, self.ny, self.nx).astype(np.float64) * scale
        self.nb, (self.cx, self.cy) = geom["cells_per_b"], geom["centre_lattice"]
        self.solid = np.zeros((self.nz, self.ny, self.nx), bool)
        for k in range(self.nz):
            self.solid[k] = k < hmap[self.y0:self.y0 + self.ny, self.x0:self.x0 + self.nx]

    def stencil(self, x_b, y_b, z_b):
        i, j, k = self.cx + self.nb * x_b - self.x0, self.cy + self.nb * y_b - self.y0, 0.5 + self.nb * z_b
        if not (0 <= i < self.nx - 1 and 0 <= j < self.ny - 1 and 0 <= k < self.nz - 1):
            raise ValueError("point outside the box: %s" % ((x_b, y_b, z_b),))
        i0, j0, k0 = int(i), int(j), int(k)
        fi, fj, fk = i - i0, j - j0, k - k0
        w = np.array([1 - fk, fk])[:, None, None] * np.array([1 - fj, fj])[None, :, None] * np.array([1 - fi, fi])[None, None, :]
        return (slice(k0, k0 + 2), slice(j0, j0 + 2), slice(i0, i0 + 2)), w

    def sample(self, x_b, y_b, z_b, field=None):
        sl, w = self.stencil(x_b, y_b, z_b)
        a = self.a if field is None else field
        return (a[(slice(None),) + sl] * w[None]).sum(axis=(1, 2, 3)), bool(self.solid[sl].any())

    def centreline_u(self, z_b, field=None):
        """Mean U along y/b = 0 at height z/b, with x/b of every cell column of the box."""
        j, k = self.cy - self.y0, 0.5 + self.nb * z_b
        j0, k0 = int(j), int(k)
        fj, fk = j - j0, k - k0
        u = (self.a if field is None else field)[0]
        col = (u[k0, j0] * (1 - fj) + u[k0, j0 + 1] * fj) * (1 - fk) + (u[k0 + 1, j0] * (1 - fj) + u[k0 + 1, j0 + 1] * fj) * fk
        return (np.arange(self.nx) + self.x0 - self.cx) / self.nb, col


def reattach_sim(box, z_b, line=None):
    xs, u = box.centreline_u(z_b) if line is None else line
    neg = (xs > 0.55) & (u < 0)
    if not neg.any():
        return None
    last = int(np.where(neg)[0].max())
    if last + 1 >= len(xs):
        return None
    return float(xs[last] + (xs[last + 1] - xs[last]) * (-u[last]) / (u[last + 1] - u[last]) - 0.5)


def reattach_tunnel(points, u_h, z_b):
    pts = sorted((p["x_b"], p["U"] / u_h) for p in points if p["section"] == "vertical" and p["z_b"] == z_b and p["x_b"] > 0.5)
    for (xa, ua), (xb, ub) in zip(pts, pts[1:]):
        if ua < 0 <= ub:
            return float(xa + (xb - xa) * (-ua) / (ub - ua) - 0.5)
    return None


def jackknife_se(values):
    """Delete-one-batch jackknife: values[k] is the statistic without batch k."""
    v = np.asarray(values, float)
    n = len(v)
    return float(np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum()))


def group_stats(ps):
    out = {"n": len(ps)}
    if not ps:
        return out
    for q in ("Uh", "U", "V", "W"):
        P, O = np.array([p["pred"][q] for p in ps]), np.array([p["obs"][q] for p in ps])
        hits = int(sum(is_hit(a, b) for a, b in zip(P, O)))
        d = P - O
        out[q] = {"hits": hits, "hit_rate": hits / len(ps), "interval95": wilson(hits, len(ps)), "bias": float(d.mean()), "rmse": float(np.sqrt((d ** 2).mean()))}
        if q in ("Uh", "U"):  # sampling error of the rate itself: the rate without one batch, for every batch
            loo = [p.get("_loo_hit_" + q) for p in ps]
            out[q]["standard_error"] = None if loo[0] is None else jackknife_se(np.mean(np.array(loo, float), axis=0))
            se = [p["standard_error"][q] for p in ps if p.get("standard_error")]
            out[q]["mean_standard_error_of_the_points"] = float(np.mean(se)) if se else None
    P, O = np.array([p["pred"]["Uh"] for p in ps]), np.array([p["obs"]["Uh"] for p in ps])
    ok = O > 0
    out["fac2_Uh"] = float(np.mean((P[ok] >= 0.5 * O[ok]) & (P[ok] <= 2.0 * O[ok]))) if ok.any() else None
    if "sigma" in ps[0]:
        S, T = np.array([p["sigma"] for p in ps]), np.array([[p["obs"][k] for k in ("su", "sv", "sw")] for p in ps])
        out["sigma"] = {"simulation_su_sv_sw": S.mean(axis=0).tolist(), "tunnel_su_sv_sw": T.mean(axis=0).tolist(), "share_of_tunnel_su_sv_sw": (S.mean(axis=0) / T.mean(axis=0)).tolist()}
    return out


def main():
    run, geom_dir, meas_path, out_path = (Path(a).resolve() for a in sys.argv[1:5])
    planes = run / "planes"
    pm = dict(ln.split() for ln in (planes / "meta.txt").read_text().splitlines())
    geom = json.loads((geom_dir / "meta.json").read_text())
    meas = json.loads(meas_path.read_text())
    u_ref, u_h = float(pm["u_ref_lattice"]), meas["U_H_ms"]
    hmap = np.fromfile(geom_dir / "height.u16", np.uint16).reshape(geom["ny"], geom["nx"])
    locs, pairs = unique_locations(meas["points"], u_h)
    box = Box(planes, "vol_mean.f32", "vol_meta.txt", 3, 1.0 / u_ref, geom, hmap)
    std = Box(planes, "vol_std.f32", "vol_meta.txt", 4, 1.0, geom, hmap) if (planes / "vol_std.f32").exists() else None
    var = None if std is None else (std.a[:3] / u_ref) ** 2  # variances in U_H^2: these are interpolated, then the root is taken

    for p in locs:
        u, touch = box.sample(p["x_b"], p["y_b"], p["z_b"])
        p["pred"] = {"U": float(u[0]), "V": float(u[1]), "W": float(u[2]), "Uh": math.hypot(u[0], u[1])}
        p["touches_solid"] = touch
        p["hit_Uh"], p["hit_U"] = is_hit(p["pred"]["Uh"], p["obs"]["Uh"]), is_hit(p["pred"]["U"], p["obs"]["U"])
        if var is not None:
            p["sigma"] = np.sqrt(np.maximum(std.sample(p["x_b"], p["y_b"], p["z_b"], var)[0], 0.0)).tolist()

    # batches of the averaging window from the means so far of the frames, and what they give: standard errors per location,
    # hit rates and reattachment lengths without one batch (jackknife)
    tagged = {}
    for f in planes.glob("vol_mean-*.f32"):
        tag = f.stem.split("-")[1]
        if (planes / ("vol_meta-%s.txt" % tag)).exists():
            b = Box(planes, f.name, "vol_meta-%s.txt" % tag, 3, 1.0 / u_ref, geom, hmap)
            if 0 < b.samples < box.samples:
                tagged[b.samples] = b
    n_all = box.samples
    cum = sorted(tagged.items()) + [(n_all, box)]
    cand = [(n - (cum[i - 1][0] if i else 0), (b.a * n - cum[i - 1][1].a * cum[i - 1][0]) if i else b.a * n) for i, (n, b) in enumerate(cum)]  # (samples, sum field) per batch
    lens = [c[0] for c in cand]
    n_batch = max(set(lens), key=lens.count)
    batches = [c[1] / c[0] for c in cand if c[0] == n_batch] if sum(c[0] == n_batch for c in cand) == len(cand) and len(cand) >= 3 else []
    nbt = len(batches)
    loo_lines = None
    if nbt:
        for p in locs:
            P = np.array([box.sample(p["x_b"], p["y_b"], p["z_b"], b)[0] for b in batches])            # [batch, component]
            uh = np.hypot(P[:, 0], P[:, 1])
            p["standard_error"] = {"U": float(P[:, 0].std(ddof=1) / math.sqrt(nbt)), "Uh": float(uh.std(ddof=1) / math.sqrt(nbt)), "W": float(P[:, 2].std(ddof=1) / math.sqrt(nbt))}
            loo = (P.sum(axis=0)[None, :] - P) / (nbt - 1)                                              # mean of all batches but one
            p["_loo_hit_Uh"] = [is_hit(math.hypot(l[0], l[1]), p["obs"]["Uh"]) for l in loo]
            p["_loo_hit_U"] = [is_hit(l[0], p["obs"]["U"]) for l in loo]
        loo_lines = {}
        for z in (0.125, 0.5, 1.0, 1.5):
            lines = np.array([box.centreline_u(z, b)[1] for b in batches])
            xs = box.centreline_u(z)[0]
            loo_lines[z] = [reattach_sim(box, z, (xs, (lines.sum(axis=0) - lines[k]) / (nbt - 1))) for k in range(nbt)]

    groups = {"all": locs}
    for r in REGIONS[:5]:
        groups[r] = [p for p in locs if p["region"] == r]
    groups["pedestrian"] = [p for p in locs if p["pedestrian"]]
    regions = {g: group_stats(ps) for g, ps in groups.items()}
    behind = {"all_behind_block": group_stats([p for p in locs if p["x_b"] > 0.5]),
              "pedestrian_sheet_behind_block": group_stats([p for p in locs if p["x_b"] > 0.5 and p["z_b"] == 0.125]),
              "z1.25_sheet_behind_block": group_stats([p for p in locs if p["x_b"] > 0.5 and p["z_b"] == 1.25])}

    reatt = {"z_b=%g" % z: {"simulation": reattach_sim(box, z), "tunnel": reattach_tunnel(meas["points"], u_h, z),
                            "standard_error": jackknife_se(loo_lines[z]) if loo_lines and all(v is not None for v in loo_lines[z]) else None}
             for z in (0.125, 0.5, 1.0, 1.5)}
    xs, uc = box.centreline_u(0.125)
    min_u = {"simulation": float(uc[xs > 0.55].min()), "simulation_x_b": float(xs[xs > 0.55][uc[xs > 0.55].argmin()]),
             "tunnel": min(p["obs"]["U"] for p in locs if p["z_b"] == 0.125 and p["y_b"] == 0 and p["x_b"] > 0.5)}

    # settling: means so far of earlier frames -> means of parts of the averaging window

    def part(n_a, n_b):
        """Mean field over the samples n_a+1..n_b from the means so far (n_a = 0: from the start of the window)."""
        hi = box if n_b == n_all else tagged[n_b]
        return hi.a if n_a == 0 else (hi.a * n_b - tagged[n_a].a * n_a) / (n_b - n_a)

    def settle(first, second):
        d_uh, d_u = [], []
        for p in locs:
            a, b = box.sample(p["x_b"], p["y_b"], p["z_b"], first)[0], box.sample(p["x_b"], p["y_b"], p["z_b"], second)[0]
            d_uh.append(abs(math.hypot(b[0], b[1]) - math.hypot(a[0], a[1])))
            d_u.append(abs(b[0] - a[0]))
        return {"mean_abs_diff_Uh": float(np.mean(d_uh)), "max_abs_diff_Uh": float(np.max(d_uh)), "mean_abs_diff_U": float(np.mean(d_u)), "max_abs_diff_U": float(np.max(d_u))}

    settling = {"samples": n_all, "means_so_far_at_samples": sorted(tagged)}
    if n_all % 2 == 0 and n_all // 2 in tagged:
        settling["halves"] = dict(settle(part(0, n_all // 2), part(n_all // 2, n_all)), samples_each=n_all // 2)
        if n_all % 4 == 0 and 3 * n_all // 4 in tagged:
            settling["quarters"] = dict(settle(part(n_all // 2, 3 * n_all // 4), part(3 * n_all // 4, n_all)), samples_each=n_all // 4)
    judged = "quarters" if "quarters" in settling else ("halves" if "halves" in settling else None)
    settling["row_B6_judged_on"] = judged
    se_uh = [p["standard_error"]["Uh"] for p in locs if p.get("standard_error")]
    settling["as_proposed_standard_error_of_the_mean"] = ({"mean_over_locations": float(np.mean(se_uh)), "largest": float(np.max(se_uh)), "proposed_mark": 0.025, "batches": nbt,
                                                            "samples_per_batch": n_batch, "verdict": "none: the replacement rule is judged in docs/VALIDATION.md"} if se_uh else None)

    def row(value, mark, ok, **more):
        return dict({"value": value, "pass_mark": mark, "pass": ok}, **more)

    a, ped, wk = regions["all"], regions["pedestrian"], regions["wake_below_roof"]
    worst = min(B4_REGIONS, key=lambda r: regions[r]["Uh"]["hit_rate"])
    rl = reatt["z_b=0.125"]["simulation"]
    rows = {
        "B1": row(a["Uh"]["hit_rate"], ">= %.2f" % MARKS["B1"], a["Uh"]["hit_rate"] >= MARKS["B1"], interval95=a["Uh"]["interval95"], standard_error=a["Uh"]["standard_error"], n=a["n"],
                  signed_U=a["U"]["hit_rate"], signed_U_interval95=a["U"]["interval95"]),
        "B2": row(ped["Uh"]["hit_rate"], ">= %.2f" % MARKS["B2"], ped["Uh"]["hit_rate"] >= MARKS["B2"], interval95=ped["Uh"]["interval95"], standard_error=ped["Uh"]["standard_error"], n=ped["n"],
                  signed_U=ped["U"]["hit_rate"], signed_U_interval95=ped["U"]["interval95"]),
        "B3": row(wk["U"]["hit_rate"], ">= %.2f" % MARKS["B3"], wk["U"]["hit_rate"] >= MARKS["B3"], interval95=wk["U"]["interval95"], standard_error=wk["U"]["standard_error"], n=wk["n"], horizontal_speed=wk["Uh"]["hit_rate"]),
        "B4": row(regions[worst]["Uh"]["hit_rate"], "no region below %.2f" % MARKS["B4"], regions[worst]["Uh"]["hit_rate"] >= MARKS["B4"], lowest_region=worst, interval95=regions[worst]["Uh"]["interval95"],
                  standard_error=regions[worst]["Uh"]["standard_error"],
                  by_region={r: regions[r]["Uh"]["hit_rate"] for r in B4_REGIONS}, regions_below=[r for r in B4_REGIONS if regions[r]["Uh"]["hit_rate"] < MARKS["B4"]]),
        "B5": row(rl, "%.2f to %.2f b" % MARKS["B5"], rl is not None and MARKS["B5"][0] <= rl <= MARKS["B5"][1], tunnel=reatt["z_b=0.125"]["tunnel"], standard_error=reatt["z_b=0.125"]["standard_error"]),
        "B6": row(None if judged is None else [settling[judged]["mean_abs_diff_Uh"], settling[judged]["max_abs_diff_Uh"]], "mean <= %.2f and largest <= %.2f" % MARKS["B6"],
                  None if judged is None else bool(settling[judged]["mean_abs_diff_Uh"] <= MARKS["B6"][0] and settling[judged]["max_abs_diff_Uh"] <= MARKS["B6"][1]), judged_on=judged,
                  verdict_scope="as written; not a pass of the run either way", as_proposed=settling["as_proposed_standard_error_of_the_mean"]),
    }
    res = {"run": run.name, "geom": geom_dir.name, "dx_m": float(pm["dx_m"]), "cells_per_b": geom["cells_per_b"], "u_ref_lattice": u_ref, "t_si": float(pm["t_si"]), "samples": n_all,
           "avg_start_s": float(pm.get("avg_start_s", "nan")), "avg_end_s": float(pm.get("avg_end_s", "nan")), "U_H_tunnel_ms": u_h, "locations": len(locs),
           "locations_touching_solid": sum(p["touches_solid"] for p in locs),
           "rule_for_locations_measured_twice": "tunnel value = mean of the two measurements", "interval": "Wilson score, 95%",
           "rows": rows, "regions": regions, "behind_block": behind, "reattachment_behind_leeward_face_b": reatt, "min_centreline_U_pedestrian": min_u,
           "settling": settling, "batches": {"count": nbt, "samples_each": n_batch if nbt else None}, "tunnel_repeatability_pairs": pairs,
           "tunnel_repeatability_largest_difference": {k: max(pr["difference"][k] for pr in pairs) for k in ("U", "V", "W", "su", "sv", "sw")} if pairs else None,
           "fluctuation_recorded": var is not None, "points": locs}
    for p in locs:                                    # working lists of the jackknife: not part of the result
        p.pop("_loo_hit_Uh", None)
        p.pop("_loo_hit_U", None)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(res, indent=1))

    def f2(v):
        return "  -  " if v is None else "%.2f" % v

    def se2(v):
        return "     " if v is None else "±%.2f" % v

    print("%s: %d cells across b, %d samples, %d locations (%d touch a solid cell); ± = standard error from %s" %
          (run.name, geom["cells_per_b"], n_all, len(locs), res["locations_touching_solid"], ("%d batches of %d samples" % (nbt, n_batch)) if nbt else "batches (none: fewer than three equal batches)"))
    print("%-18s %4s  %-28s %-28s %6s %7s" % ("region", "n", "q(Uh) [95% interval] ±se", "q(U) [95% interval] ±se", "FAC2", "bias Uh"))
    for g in ("all",) + REGIONS:
        s = regions[g]
        if s["n"]:
            print("%-18s %4d  %.2f [%.2f, %.2f] %s      %.2f [%.2f, %.2f] %s      %6s %+7.3f" % (g, s["n"], s["Uh"]["hit_rate"], *s["Uh"]["interval95"], se2(s["Uh"]["standard_error"]),
                                                                                         s["U"]["hit_rate"], *s["U"]["interval95"], se2(s["U"]["standard_error"]), f2(s["fac2_Uh"]), s["Uh"]["bias"]))
    print("reattachment behind the leeward face [b]: " + ", ".join("%s: %s %s (tunnel %s)" % (k, f2(v["simulation"]), se2(v["standard_error"]).strip(), f2(v["tunnel"])) for k, v in reatt.items()))
    for k in ("quarters", "halves"):
        if k in settling:
            print("settling, %s of the averaging window (%d samples each): mean |dUh| %.4f, largest %.4f" % (k, settling[k]["samples_each"], settling[k]["mean_abs_diff_Uh"], settling[k]["max_abs_diff_Uh"]))
    ap = settling["as_proposed_standard_error_of_the_mean"]
    if ap:
        print("row B6 as proposed (no verdict): standard error of the mean horizontal speed, average over the locations %.4f, largest %.4f (proposed mark 0.025)" % (ap["mean_over_locations"], ap["largest"]))
    if var is not None:
        for name, s in (("all locations", regions["all"]), ("behind the block", behind["all_behind_block"]), ("  at z/b 0.125", behind["pedestrian_sheet_behind_block"]), ("  at z/b 1.25", behind["z1.25_sheet_behind_block"])):
            g = s["sigma"]
            print("fluctuation, %-17s sigma_u %.3f (tunnel %.3f, share %.2f)  sigma_v %.3f (%.3f)  sigma_w %.3f (%.3f)" %
                  (name + ":", g["simulation_su_sv_sw"][0], g["tunnel_su_sv_sw"][0], g["share_of_tunnel_su_sv_sw"][0], g["simulation_su_sv_sw"][1], g["tunnel_su_sv_sw"][1], g["simulation_su_sv_sw"][2], g["tunnel_su_sv_sw"][2]))
    for k, r in rows.items():
        v = r["value"]
        txt = "-" if v is None else (", ".join("%.3f" % x for x in v) if isinstance(v, list) else "%.3f" % v)
        tail = "no verdict" if r["pass"] is None else ("PASS" if r["pass"] else "FAIL")
        if k == "B6" and r["pass"] is not None:
            tail = "met as written" if r["pass"] else "NOT MET AS WRITTEN"
        if r.get("standard_error") is not None:
            txt += " ±%.3f" % r["standard_error"]
        print("row %s: %s  (pass mark %s)  %s" % (k, txt, r["pass_mark"], tail))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
