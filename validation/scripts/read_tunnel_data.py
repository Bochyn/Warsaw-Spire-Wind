"""Parse the AIJ benchmark case A workbook (CaseA_1_1_2.xls) into JSON and fit the approach-flow exponent.

usage: uv run --with xlrd python -I read_tunnel_data.py WORKBOOK_XLS OUT_JSON
The workbook is untrusted input (downloaded file): this script only reads cells with xlrd; it never executes anything.
Sheets (as found on 2026-11-09, printed with inspect, not assumed):
  'Geometry&Points'                 figure captions only (embedded pictures, no numbers)
  'Inflow'                          z/b, U(m/s), sigma_u, sigma_v, sigma_w, k at 24 heights 0.0625b..7.5b
  'Results (Ver. sec.)'             66 points in the vertical centre plane y/b = 0; column 'z ' is z/b
  'Results(Hol. sec.; z=0.125b)'    60 points at z/b = 0.125 (pedestrian level), y/b 0..-2 (one side)
  'Results (Hol. sec.; z=1.25b)'    60 points at z/b = 1.25; the z/b column wrongly repeats 0.125, the sheet name
                                    and the different values show it is a different height; we use 1.25.
All speeds are m/s in the tunnel; the reference is the inflow speed at building height z = 2b.
"""
import json
import math
import sys

import xlrd

SRC, OUT = sys.argv[1], sys.argv[2]
wb = xlrd.open_workbook(SRC)
sheets = {s.name.replace("）", ")"): s for s in wb.sheets()}  # the two horizontal sheets end with a full-width ')'


def rows(sheet, ncol):
    """Numeric rows of a sheet: the first ncol+1 columns (column 0 is the point number, '' in Inflow)."""
    out = []
    for r in range(sheet.nrows):
        v = sheet.row_values(r)[: ncol + 1]
        if len(v) == ncol + 1 and all(isinstance(x, float) for x in v[1:]):
            out.append(v)
    return out


inflow = [{"z_b": r[1], "U": r[2], "su": r[3], "sv": r[4], "sw": r[5], "k": r[6]} for r in rows(sheets["Inflow"], 6)]
U_H = next(p["U"] for p in inflow if abs(p["z_b"] - 2.0) < 1e-9)

points = []
for sheet, section, z_fixed in (("Results (Ver. sec.)", "vertical", None), ("Results(Hol. sec.; z=0.125b)", "horizontal_0.125b", 0.125),
                                ("Results (Hol. sec.; z=1.25b)", "horizontal_1.25b", 1.25)):
    for r in rows(sheets[sheet], 10):
        z_b = r[3] if z_fixed is None else z_fixed
        points.append({"sheet": sheet, "section": section, "n": int(r[0]), "x_b": r[1], "y_b": r[2], "z_b": z_b, "z_b_in_file": r[3],
                       "U": r[4], "V": r[5], "W": r[6], "su": r[7], "sv": r[8], "sw": r[9], "k": r[10]})

# power-law fit U/U_H = (z/H)^alpha, H = 2b, in log space; amplitude pinned at U_H (the solver pins the inflow at z_ref)
def fit_alpha(zmin_b, zmax_b):
    xs, ys = [], []
    for p in inflow:
        if zmin_b <= p["z_b"] <= zmax_b:
            xs.append(math.log(p["z_b"] / 2.0))
            ys.append(math.log(p["U"] / U_H))
    a_pinned = sum(x * y for x, y in zip(xs, ys)) / sum(x * x for x in xs)
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    a_free = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    amp_free = math.exp(my - a_free * mx)
    rms = math.sqrt(sum((math.exp(a_pinned * x) - math.exp(y)) ** 2 for x, y in zip(xs, ys)) / n)
    return {"z_range_b": [zmin_b, zmax_b], "n": n, "alpha_pinned": a_pinned, "alpha_free": a_free, "amplitude_free_over_UH": amp_free, "rms_pinned": rms}


fits = [fit_alpha(0.0625, 7.5), fit_alpha(0.125, 3.5), fit_alpha(0.0625, 2.0), fit_alpha(0.125, 5.0)]
dups = {}
for p in points:
    dups.setdefault((p["x_b"], p["y_b"], p["z_b"]), []).append(p)
duplicates = [[(q["section"], q["U"], q["V"], q["W"]) for q in v] for v in dups.values() if len(v) > 1]

out = {"source": SRC.split("/")[-1], "b_m_tunnel": 0.08, "H_over_b": 2.0, "U_H_ms": U_H, "inflow": inflow, "points": points,
       "alpha_fits": fits, "duplicate_locations": len(duplicates), "unique_locations": len(dups)}
with open(OUT, "w") as f:
    json.dump(out, f, indent=1)

print("sheets:", list(sheets))
print("inflow: %d heights, z/b %g..%g, U_H = U(z/b=2) = %.3f m/s, Iu at H = %.1f%%" % (len(inflow), inflow[0]["z_b"], inflow[-1]["z_b"], U_H, 100 * next(p["su"] for p in inflow if p["z_b"] == 2.0) / U_H))
for s in ("vertical", "horizontal_0.125b", "horizontal_1.25b"):
    ps = [p for p in points if p["section"] == s]
    print("%-18s %d points, x/b %s, y/b %s, z/b %s" % (s, len(ps), sorted({p["x_b"] for p in ps}), sorted({p["y_b"] for p in ps}), sorted({p["z_b"] for p in ps})))
print("points %d, unique locations %d, duplicates (vertical vs horizontal 0.125b):" % (len(points), len(dups)))
for d in duplicates:
    print("  ", d)
for f in fits:
    print("alpha fit z/b %s: pinned %.3f (rms %.3f U_H), free %.3f with amplitude %.3f U_H" % (f["z_range_b"], f["alpha_pinned"], f["rms_pinned"], f["alpha_free"], f["amplitude_free_over_UH"]))
print("wrote", OUT)
