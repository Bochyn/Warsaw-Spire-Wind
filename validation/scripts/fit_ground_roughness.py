"""Ground roughness length for the wall-stress ground: log-law fit to the lowest heights of the tunnel's inflow sheet.

usage: uv run python -I fit_ground_roughness.py TUNNEL_JSON OUT_JSON
Rule (written down before any run with the wall-stress ground): least squares of U(z) = (u_tau/kappa) ln(z/z0),
kappa = 0.4, on the lowest THREE heights of the sheet (z/b = 0.0625, 0.125, 0.25: 2, 4, 8 m at the analogue scale b = 32 m).
Only the approach-flow sheet is read; the measuring points around the block are not touched.
The fits on more heights and the slopes between neighbouring heights are printed as a sensitivity, not used.
"""
import json
import math
import sys

MEAS, OUT = sys.argv[1], sys.argv[2]
KAPPA, B_M = 0.4, 32.0
m = json.load(open(MEAS))
UH = m["U_H_ms"]
rows = sorted(m["inflow"], key=lambda p: p["z_b"])


def fit(n):
    xs = [math.log(p["z_b"] * B_M) for p in rows[:n]]
    ys = [p["U"] / UH for p in rows[:n]]
    mx, my = sum(xs) / n, sum(ys) / n
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    z0 = math.exp(mx - my / slope)  # U = slope * (ln z - ln z0)
    res = [slope * (x - math.log(z0)) - y for x, y in zip(xs, ys)]
    return {"heights_z_b": [p["z_b"] for p in rows[:n]], "u_tau_over_UH": KAPPA * slope, "z0_m": z0, "max_abs_residual_UH": max(abs(r) for r in res)}


def cd(z0, z1):
    return (KAPPA / math.log(z1 / z0)) ** 2


rule = fit(3)
out = {"rule": "least squares of U = (u_tau/kappa) ln(z/z0), kappa = 0.4, lowest three heights of the Inflow sheet", "kappa": KAPPA, "b_m": B_M,
       "z0_m": rule["z0_m"], "u_tau_over_UH": rule["u_tau_over_UH"], "fit": rule,
       "Cd_first_cell": {"16 cells per b (z1 = 1 m)": cd(rule["z0_m"], 1.0), "24 cells per b (z1 = 0.667 m)": cd(rule["z0_m"], 2.0 / 3.0), "32 cells per b (z1 = 0.5 m)": cd(rule["z0_m"], 0.5)},
       "sensitivity_not_used": {"fit_on_n_lowest_heights": {str(n): fit(n) for n in (2, 4, 5, 6)},
                                "u_tau_over_UH_from_neighbouring_heights": [
                                    {"z_b": [a["z_b"], b["z_b"]], "u_tau_over_UH": KAPPA * (b["U"] - a["U"]) / UH / math.log(b["z_b"] / a["z_b"])} for a, b in zip(rows[:5], rows[1:6])],
                                "u_tau_over_UH_if_sigma_u_were_2.45_u_tau": [{"z_b": p["z_b"], "u_tau_over_UH": p["su"] / 2.45 / UH} for p in rows[:4]]}}
json.dump(out, open(OUT, "w"), indent=1)
print("rule (3 lowest heights): z0 = %.3e m, u_tau = %.4f U_H, max residual %.4f U_H" % (rule["z0_m"], rule["u_tau_over_UH"], rule["max_abs_residual_UH"]))
for k, v in out["Cd_first_cell"].items():
    print("  Cd, %s: %.5f" % (k, v))
for n, f in out["sensitivity_not_used"]["fit_on_n_lowest_heights"].items():
    print("  not used, %s heights: z0 = %.3e m, u_tau = %.4f U_H, Cd(1 m) = %.5f" % (n, f["z0_m"], f["u_tau_over_UH"], cd(f["z0_m"], 1.0)))
print("  not used, u_tau/U_H from neighbouring heights:", ["%.4f" % d["u_tau_over_UH"] for d in out["sensitivity_not_used"]["u_tau_over_UH_from_neighbouring_heights"]])
print("  not used, u_tau/U_H if sigma_u = 2.45 u_tau:", ["%.4f" % d["u_tau_over_UH"] for d in out["sensitivity_not_used"]["u_tau_over_UH_if_sigma_u_were_2.45_u_tau"]])
