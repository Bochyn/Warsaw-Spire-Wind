"""Inflow table with a calibrated length scale and, optionally, a correlation of u and w at the inlet. A wrapper:
make_inflow_table.py stays the tool for gains alone.

usage: uv run --with numpy python -I make_calibrated_inflow.py CELLS_PER_B OUT_DIR TUNNEL_JSON GAIN_JSON|- L_FACTOR [RHO_JSON]
 TUNNEL_JSON the output of read_tunnel_data.py
 GAIN_JSON  the calibration file ({"z_b": [...], "su": [...], "sv": [...], "sw": [...]}, optional "U"), "-" for none
 L_FACTOR   the streamwise integral scale is L_FACTOR times the length of the rule in make_inflow_table.py
            (1.0 as declared; 2.6 is the integral scale consistent with the same k and epsilon)
 RHO_JSON   optional: {"z_b": [...], "rho": [...]}, the correlation coefficient of u and w by height (negative for a
            downward momentum flux), interpolated linearly in z/b and held outside the range (see correlate_uw).
            Not used for the tables of the wind-tunnel check.
Writes OUT_DIR/inflow_b<NB>.bin and .json exactly as make_inflow_table.py does (the sidecar records L_factor and the correlation).
Importing this file runs nothing: correlate_uw() can be used by other cases.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inflow_table as mi  # noqa: E402


def correlate_uw(tab, rho, su, sw):
    """Optional correlation of u and w: the table's w becomes sqrt(1 - rho^2) w + rho (sigma_w / sigma_u) u per layer, so sigma_w stays and the
    covariance of u and w becomes rho sigma_u sigma_w. tab: a table as make_inflow_table.build() or read_table() returns it, changed
    in place; rho, su, sw: per layer (the correlation coefficient, limited to +-0.95, and the two sigmas the table imposes).
    The playback formula is linear in the stored factors B, so this needs no change of the file format or the solver."""
    rho = np.clip(rho, -0.95, 0.95)
    mix = np.where(su > 0, rho * sw / np.maximum(su, 1e-12), 0.0)  # share of u added to w, per layer
    keep = np.sqrt(1.0 - rho ** 2)
    for key in ("B_re", "B_im"):
        b = tab[key].astype(np.float64)
        b[2] = keep[:, None] * b[2] + mix[:, None] * b[0]
        tab[key] = b.astype(np.float32)


def main():
    nb, out, sheet, gain_arg, l_factor = int(sys.argv[1]), Path(sys.argv[2]).resolve(), sys.argv[3], sys.argv[4], float(sys.argv[5])
    rho_path = Path(sys.argv[6]) if len(sys.argv) > 6 else None
    if not 1.0 <= l_factor <= 2.6:
        raise SystemExit("L_FACTOR must lie between 1.0 and 2.6")
    mi.L_FACTOR = l_factor  # read by profiles() and build() when they run
    sys.argv = [sys.argv[0], str(nb), str(out), sheet] + ([gain_arg] if gain_arg != "-" else [])
    mi.main()
    if rho_path is not None:
        rho_in = json.loads(rho_path.read_text())
        path, side = out / ("inflow_b%d.bin" % nb), out / ("inflow_b%d.json" % nb)
        t, meta = mi.read_table(path), json.loads(side.read_text())
        z = np.array(meta["layers"]["z_b"])
        su, sw = np.array(meta["layers"]["su"]), np.array(meta["layers"]["sw"])
        if meta.get("gain"):
            su, sw = su * np.interp(z, meta["gain"]["z_b"], meta["gain"]["su"]), sw * np.interp(z, meta["gain"]["z_b"], meta["gain"]["sw"])
        rho = np.interp(z, rho_in["z_b"], rho_in["rho"])
        correlate_uw(t, rho, su, sw)
        mi.write_table(path, t)
        meta["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        meta["uw_correlation"] = {"z_b": rho_in["z_b"], "rho": rho_in["rho"], "rule": "w' = sqrt(1 - rho^2) w + rho (sigma_w / sigma_u) u per layer"}
        side.write_text(json.dumps(meta, indent=1))
        print("u-w correlation written into the table: rho from %.2f to %.2f; new sha256 %s" % (rho.min(), rho.max(), meta["sha256"][:16]))


if __name__ == "__main__":
    main()
