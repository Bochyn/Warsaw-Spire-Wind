"""Height map and grid for the AIJ case A analogue at full scale: one block b x b x 2b on flat ground, wind along +x.

usage: uv run --with numpy python -I make_block_geometry.py CELLS_PER_B OUT_DIR
Full-scale analogue: b = 32 m, H = 2b = 64 m (the tunnel model is 0.08 m; the geometry is scaled 400:1, so x/b, y/b, z/b
are the same). Domain (AIJ guideline, Tominaga et al. 2008): 5H upstream of the building, 10H downstream, 5H to each side
plus half a building, about 5.9H high (the top row is the boundary). Height map: little-endian uint16, HY rows of HX
columns, 1 = ground, 1 + H/dx inside the block (solid from the ground layer z = 0 up to cell H/dx; the roof face is at
height H above the ground, which sits at z = 0.5 cells). Writes height.u16 and meta.json (cell indices of the block,
its centre in continuous lattice coordinates, the 3D output box) and prints the environment line for FluidX3D.
"""
import json
import sys
from pathlib import Path

import numpy as np

NB = int(sys.argv[1])  # cells across the building width b
OUT = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
B_M, H_M = 32.0, 64.0
dx = B_M / NB
NH = 2 * NB  # cells of building height
up, down, side = 10 * NB, 20 * NB, 10 * NB + NB // 2  # 5H, 10H, 5H + b/2
HX, HY = up + NB + down, 2 * side + NB
NZ = 1 + NH + 10 * NB - NB // 2  # ground layer + building + 4.75H of air above the roof, top row = boundary
H = np.ones((HY, HX), np.uint16)
x0, y0 = up, side
H[y0:y0 + NB, x0:x0 + NB] = 1 + NH
H.tofile(OUT / "height.u16")
cx, cy = x0 + NB / 2 - 0.5, y0 + NB / 2 - 0.5  # building centre, continuous lattice index (cell centres at integers)
# 3D output box: measurement points span x/b -0.75..3.25, y/b -2..0, z/b <= 3.5; box x/b -4..5, y/b -3..3, z/b 0..4.5
box = {"x0": int(cx - 4 * NB) - 1, "x1": int(cx + 5 * NB) + 2, "y0": int(cy - 3 * NB) - 1, "y1": int(cy + 3 * NB) + 2, "z1": int(0.5 + 4.5 * NB) + 2}
meta = {"b_m": B_M, "H_m": H_M, "dx_m": dx, "cells_per_b": NB, "nx": HX, "ny": HY, "nz": NZ, "block_x": [x0, x0 + NB - 1], "block_y": [y0, y0 + NB - 1],
        "block_top_cell": NH, "centre_lattice": [cx, cy], "ground_lattice_z": 0.5, "box": box,
        "upstream_H": up / NH, "downstream_H": down / NH, "side_H": side / NH, "air_height_H": (NZ - 2) / NH,
        "blockage": (NB * NH) / (HY * (NZ - 2)), "cells_M": HX * HY * NZ / 1e6,
        "note": "point (x/b, y/b, z/b) -> lattice (cx + NB*x/b, cy + NB*y/b, 0.5 + NB*z/b); box indices are cells"}
(OUT / "meta.json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta, indent=1))
print("ENV: HX=%d HY=%d NZ=%d DX_M=%.6f PX=%d PY=%d VX0=%d VX1=%d VY0=%d VY1=%d VZ1=%d ZPED=%d ZHI=%d" %
      (HX, HY, NZ, dx, round(cx), round(cy), box["x0"], box["x1"], box["y0"], box["y1"], box["z1"], max(int(0.5 + 0.125 * NB), 1), int(0.5 + 1.25 * NB)))
