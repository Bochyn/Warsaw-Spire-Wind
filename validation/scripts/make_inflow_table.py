"""Synthetic inflow turbulence for the AIJ case A analogue: mode table and reference evaluator.

usage: uv run --with numpy python -I make_inflow_table.py CELLS_PER_B OUT_DIR [TUNNEL_JSON] [GAIN_JSON]
TUNNEL_JSON is the output of read_tunnel_data.py (default ../results/tunnel-data.json).
GAIN_JSON (optional, a calibration file such as ../inflow-calibration.json): {"z_b": [...], "su": [...], "sv": [...], "sw": [...]},
multipliers of the three sigmas by height, interpolated linearly and held outside their range; lengths and time scales
stay those of the sheet. Optional key "U" (the mean-speed correction, added 2026-11-10): multipliers of the mean speed
by height, same interpolation; without the key the mean speed is the sheet's and the table is byte for byte what it was.
The corrected mean speed is the table's U: the solver puts it on every fixed face (inlet, outlet, sides, top) and into the
initial field, so the inlet plane and the faces that meet it carry one profile and no step runs along their common edge.
Everything else in the table (modes, phases, amplitudes, the time scale per height) is built from the sheet's own speed
and does not change with "U".
Writes OUT_DIR/inflow_b<NB>.bin (the table the solver plays back) and OUT_DIR/inflow_b<NB>.json (profiles, derived
length scales, every parameter, sha256 of the table) and prints a summary line for the run record.
File format: write_table() below; the solver reads it with read_inflow() in src/setup.cpp. Everything is in units of the building
width b and of the approach speed at building height U_H (time in b/U_H); the full-scale analogue is b = 32 m,
U_H = 10 m/s, the tunnel b = 0.08 m, U_H = 4.491 m/s.

The fluctuation on the inlet plane x = 0 is a sum of N random Fourier modes
    u_i(y, z, t) = wy(y) wz(z) g_i(z) sum_n q_n(z) s_ni cos(ky_n y + kz_n z + psi_n - Om_n t)   (i = u, v, w)
followed by the removal of the taper-weighted plane mean of u. The table stores the factors of this formula:
    B_i[z, n] = wz(z) g_i(z) q_n(z) s_ni exp(i kz_n z),   Ey[n, y] = exp(i ky_n y),   Om_n, psi_n, wy, wz, U(z)
so that u_i = wy(y) Re( sum_n B_i[z, n] Ey[n, y] exp(i (psi_n - Om_n t)) ).
"""
from __future__ import annotations

import hashlib
import json
import struct
import sys
from pathlib import Path

import numpy as np

# ---- declared parameters ----
SEED = 20261110            # date of the declaration; not tuned
N_MASTER = 512             # modes drawn for the finest declared grid; a coarser grid keeps those it can carry
NB_MASTER = 32             # finest declared grid, cells across b
LAMBDA_MIN_CELLS = 4.0     # shortest wavelength in cells: k_max = 2 pi / (4 dx)
K_LO_OVER_KE = 0.1         # lowest wavenumber as a fraction of the reference energy wavenumber
C_MU = 0.09                # constant of the length-scale rule (sheet_length)
FLOOR_CELLS = 8.0          # floor of the length scale in cells
L_FACTOR = 1.0             # streamwise integral scale = L_FACTOR * length of the rule (switch; 2.6 would be Pope's L11 for the same k, eps)
SLOPE_HALF_WINDOW = 2      # dU/dz from a least-squares line of U against ln z over 2*2+1 sheet heights
U_C = 1.0                  # convection speed of the pattern in U_H (the reference speed)
TAPER_CELLS = 8            # fluctuation rises from 0 at a side or top face cell to full within 8 cells
MAGIC = b"SCINFLW1"


def load_sheet(path):
    d = json.loads(Path(path).read_text())
    uh = float(d["U_H_ms"])
    inf = d["inflow"]
    col = lambda key: np.array([float(r[key]) for r in inf])
    return {"z": col("z_b"), "U": col("U") / uh, "su": col("su") / uh, "sv": col("sv") / uh, "sw": col("sw") / uh,
            "k": col("k") / uh ** 2, "U_H_ms": uh}


def grid(nb):
    """Inlet plane of make_geom.py: ny cells across, nz layers (0 = ground, nz-1 = top face)."""
    side = 10 * nb + nb // 2
    ny, nz = 2 * side + nb, 1 + 2 * nb + 10 * nb - nb // 2
    y = (np.arange(ny) - (ny - 1) / 2.0) / nb          # centre of cell iy, in b, origin on the plane's centre line
    z = (np.arange(nz) - 0.5) / nb                      # centre of layer j above the ground (wall at lattice z = 0.5), in b
    return ny, nz, y, z


def sheet_length(sh):
    """Length-scale rule at the 24 sheet heights: L = C_mu^(3/4) k^(3/2) / eps, eps = C_mu^(1/2) k dU/dz  =>  L = C_mu^(1/4) k^(1/2) / (dU/dz)."""
    z, U = sh["z"], sh["U"]
    lz = np.log(z)
    dudz = np.empty_like(z)
    for i in range(len(z)):
        lo, hi = max(0, i - SLOPE_HALF_WINDOW), min(len(z), i + SLOPE_HALF_WINDOW + 1)
        dudz[i] = np.polyfit(lz[lo:hi], U[lo:hi], 1)[0] / z[i]
    return C_MU ** 0.25 * np.sqrt(sh["k"]) / dudz, dudz


def profiles(sh, nb):
    """Targets per cell layer. U: linear in ln z between sheet heights; sigmas and L: linear in z. Below the lowest
    sheet height the lowest segment is continued (for U a log law through the two lowest points; for the sigmas a
    straight line, so that reading the layers back at the lowest sheet height returns the sheet value). Above the
    highest sheet height everything is held at its last value."""
    ny, nz, y, z = grid(nb)
    zs = sh["z"]
    zc = np.maximum(z, 1e-9)
    U = np.interp(np.log(zc), np.log(zs), sh["U"])
    low = zc < zs[0]
    slope0 = (sh["U"][1] - sh["U"][0]) / np.log(zs[1] / zs[0])
    U[low] = np.maximum(sh["U"][0] + slope0 * np.log(zc[low] / zs[0]), 0.0)
    L_sheet, dudz = sheet_length(sh)
    out = {"z": z, "U": U, "L_rule": np.interp(z, zs, L_sheet)}
    for key in ("su", "sv", "sw"):
        out[key] = np.interp(z, zs, sh[key])
        out[key][low] = np.maximum(sh[key][0] + (sh[key][1] - sh[key][0]) / (zs[1] - zs[0]) * (zc[low] - zs[0]), 0.0)
    floor = FLOOR_CELLS / nb
    out["L"] = np.maximum(L_FACTOR * out["L_rule"], floor)   # streamwise integral scale wanted in the domain
    out["floor_binds"] = L_FACTOR * out["L_rule"] < floor
    # the pattern passes the inlet at U_C at every height, the flow then carries it at U(z): a probe sees the time
    # scale L_spec/U_C, the domain the length U(z) L_spec/U_C. L_spec makes that length equal to L.
    out["L_spec"] = out["L"] * U_C / np.maximum(U, 1e-6)
    for key in ("U", "su", "sv", "sw"):
        out[key][0] = 0.0                                    # ground layer: solid
    return out, L_sheet, dudz


def vk_shape(k, ke):
    """von Karman energy spectrum, shape only: E(k) ~ (k/ke)^4 / (1 + (k/ke)^2)^(17/6)."""
    x = k / ke
    return x ** 4 / (1.0 + x * x) ** (17.0 / 6.0)


def _trapz(f, x):
    return np.sum(0.5 * (f[..., 1:] + f[..., :-1]) * np.diff(x), axis=-1)


def l11_truncated(ke, kq):
    """Longitudinal integral scale of an isotropic field whose spectrum is the von Karman shape cut to the band kq:
    L11 = (3 pi / 4) int E/k dk / int E dk."""
    e = vk_shape(kq[None, :], np.atleast_1d(ke)[:, None])
    return 0.75 * np.pi * _trapz(e / kq, kq) / _trapz(e, kq)


def solve_ke(L_spec, kq):
    """ke(z) so that the cut spectrum has the wanted integral scale (bisection in ln ke; L11 falls as ke rises)."""
    lo = np.full(L_spec.shape, np.log(1e-3))
    hi = np.full(L_spec.shape, np.log(kq[-1]))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        too_long = l11_truncated(np.exp(mid), kq) > L_spec
        lo, hi = np.where(too_long, mid, lo), np.where(too_long, hi, mid)
    return np.exp(0.5 * (lo + hi))


def master_modes(ke_ref, us):
    """N_MASTER modes for the finest declared grid, sorted by |k|. Only raw uniform doubles of PCG64 are used, so the
    same seed gives the same modes in every numpy version."""
    r = np.random.Generator(np.random.PCG64(SEED)).random((N_MASTER, 7))
    k_lo, k_hi = K_LO_OVER_KE * ke_ref, 2.0 * np.pi / (LAMBDA_MIN_CELLS / NB_MASTER)
    kq = np.geomspace(k_lo, k_hi, 20001)
    e = vk_shape(kq, ke_ref)
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(kq))))
    cdf /= cdf[-1]
    kmag = np.interp((np.arange(N_MASTER) + r[:, 0]) / N_MASTER, cdf, kq)   # one mode per slice of equal energy
    mu = 2.0 * r[:, 1] - 1.0                                               # direction uniform on the sphere
    az = 2.0 * np.pi * r[:, 2]
    st = np.sqrt(1.0 - mu * mu)
    khat = np.stack([mu, st * np.cos(az), st * np.sin(az)], axis=1)
    e1 = np.stack([-st, mu * np.cos(az), mu * np.sin(az)], axis=1)         # unit vectors normal to khat
    e2 = np.stack([np.zeros(N_MASTER), -np.sin(az), np.cos(az)], axis=1)
    beta = 2.0 * np.pi * r[:, 3]
    s = np.cos(beta)[:, None] * e1 + np.sin(beta)[:, None] * e2            # amplitude direction, s . k = 0
    psi = 2.0 * np.pi * r[:, 4]
    xi = np.sqrt(-2.0 * np.log(1.0 - r[:, 5])) * np.cos(2.0 * np.pi * r[:, 6])   # Box-Muller, standard normal
    kvec = kmag[:, None] * khat
    omega = U_C * kvec[:, 0] + us * kmag * xi                              # convection plus random sweeping
    return {"k": kvec, "kmag": kmag, "s": s, "psi": psi, "omega": omega, "k_lo": k_lo, "k_hi_master": k_hi}


def build(sheet_path, nb, gain=None):
    sh = load_sheet(sheet_path)
    ny, nz, y, z = grid(nb)
    if nb > NB_MASTER:
        raise SystemExit("grid finer than the declared master grid: raise NB_MASTER (this changes the table)")
    prof, L_sheet, dudz = profiles(sh, nb)
    iH = int(np.argmin(np.abs(sh["z"] - 2.0)))
    ke_ref = 0.7468343 / (L_FACTOR * L_sheet[iH])                          # untruncated von Karman: ke = 0.7468 / L11
    us = float(np.sqrt(2.0 * sh["k"][iH] / 3.0))                           # one-component rms at building height
    m = master_modes(ke_ref, us)
    k_hi = 2.0 * np.pi / (LAMBDA_MIN_CELLS / nb)
    keep = m["kmag"] <= k_hi
    k, kmag, s, psi, omega = m["k"][keep], m["kmag"][keep], m["s"][keep], m["psi"][keep], m["omega"][keep]
    n = int(keep.sum())
    kq = np.geomspace(m["k_lo"], k_hi, 4001)
    ke = solve_ke(prof["L_spec"], kq)                                      # per layer
    q = np.sqrt(vk_shape(kmag[None, :], ke[:, None]) / vk_shape(kmag[None, :], ke_ref))   # [nz, n], importance weights
    raw = np.sqrt(0.5 * np.einsum("zn,ni->zi", q * q, s * s))              # exact rms of the unscaled field per component
    target = np.stack([prof["su"], prof["sv"], prof["sw"]], axis=1)
    if gain:
        target = target * np.stack([np.interp(z, gain["z_b"], gain[key]) for key in ("su", "sv", "sw")], axis=1)
    g = target / raw                                                       # [nz, 3]
    if gain and "U" in gain:
        gu = np.asarray(gain["U"], float)
        if gu.shape != np.asarray(gain["z_b"], float).shape or not np.all(np.isfinite(gu)) or not np.all(gu > 0.0):
            raise SystemExit('gain key "U": one finite multiplier above zero per height of z_b')
    u_gain = np.interp(z, gain["z_b"], gain["U"]) if gain and "U" in gain else None   # mean-speed correction: mean speed imposed by the table
    u_table = prof["U"] * u_gain if u_gain is not None else prof["U"]
    d_top = (nz - 1) - np.arange(nz)
    wz = np.where(d_top < TAPER_CELLS, np.sin(0.5 * np.pi * d_top / TAPER_CELLS) ** 2, 1.0)
    wz[0] = 0.0
    d_side = np.minimum(np.arange(ny), (ny - 1) - np.arange(ny))
    wy = np.where(d_side < TAPER_CELLS, np.sin(0.5 * np.pi * d_side / TAPER_CELLS) ** 2, 1.0)
    amp = (wz[:, None] * q)[None, :, :] * g.T[:, :, None] * s.T[:, None, :]            # [3, nz, n]
    B = amp * np.exp(1j * z[None, :, None] * k[None, None, :, 2])
    Ey = np.exp(1j * k[:, 1][:, None] * y[None, :])                                      # [n, ny]
    tab = {"nb": nb, "ny": ny, "nz": nz, "n": n, "taper": TAPER_CELLS, "seed": SEED, "u_c": U_C,
           "omega": omega.astype(np.float64), "psi": psi.astype(np.float64),
           "U": u_table.astype(np.float32), "wy": wy.astype(np.float32), "wz": wz.astype(np.float32),
           "B_re": B.real.astype(np.float32), "B_im": B.imag.astype(np.float32),
           "Ey_re": Ey.real.astype(np.float32), "Ey_im": Ey.imag.astype(np.float32)}
    design = {"sheet": sh, "prof": prof, "L_sheet": L_sheet, "dudz_sheet": dudz, "ke": ke, "ke_ref": ke_ref, "us": us,
              "k": k, "kmag": kmag, "s": s, "q": q, "g": g, "k_lo": m["k_lo"], "k_hi": k_hi, "y": y, "z": z, "u_gain": u_gain, "u_table": u_table,
              "n_eff": (q ** 2).sum(1) ** 2 / np.maximum((q ** 4).sum(1), 1e-300)}
    return tab, design


# ---- file format (little-endian; the layout is the order of the writes below) ----
def write_table(path, t):
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<6i", 1, t["nb"], t["ny"], t["nz"], t["n"], t["taper"]))
        f.write(struct.pack("<q", t["seed"]))
        f.write(struct.pack("<d", t["u_c"]))
        for key in ("omega", "psi"):
            f.write(t[key].astype("<f8").tobytes())
        for key in ("U", "wy", "wz", "B_re", "B_im", "Ey_re", "Ey_im"):
            f.write(np.ascontiguousarray(t[key]).astype("<f4").tobytes())


def read_table(path):
    raw = Path(path).read_bytes()
    if raw[:8] != MAGIC:
        raise ValueError("not an inflow table")
    ver, nb, ny, nz, n, taper = struct.unpack_from("<6i", raw, 8)
    seed, = struct.unpack_from("<q", raw, 32)
    u_c, = struct.unpack_from("<d", raw, 40)
    t = {"nb": nb, "ny": ny, "nz": nz, "n": n, "taper": taper, "seed": seed, "u_c": u_c}
    off = 48
    for key in ("omega", "psi"):
        t[key] = np.frombuffer(raw, "<f8", n, off).copy()
        off += 8 * n
    for key, shape in (("U", (nz,)), ("wy", (ny,)), ("wz", (nz,)), ("B_re", (3, nz, n)), ("B_im", (3, nz, n)), ("Ey_re", (n, ny)), ("Ey_im", (n, ny))):
        cnt = int(np.prod(shape))
        t[key] = np.frombuffer(raw, "<f4", cnt, off).reshape(shape).copy()
        off += 4 * cnt
    if ver != 1 or off != len(raw):
        raise ValueError("inflow table: version or size mismatch")
    return t


# ---- reference evaluator ----
def time_of_step(step, u_ref_lattice, nb):
    """Solver step -> time in b/U_H. The solver's dt = dx u_ref / u_si, so step dt U_H / b = step u_ref / nb."""
    return step * u_ref_lattice / nb


def plane(t, time, dtype=np.float32):
    """u'(y, z, time) on the whole inlet plane, shape [3, nz, ny], in U_H, exactly as the solver must compute it:
    phases in double, the sums in single precision, then the mean of u over the plane removed with the taper as weight."""
    ph = t["psi"] - t["omega"] * float(time)
    c, s = np.cos(ph).astype(dtype), np.sin(ph).astype(dtype)
    br, bi = t["B_re"].astype(dtype), t["B_im"].astype(dtype)
    cr, ci = br * c - bi * s, br * s + bi * c                              # B exp(i phase)
    with np.errstate(all="ignore"):                                        # Accelerate raises spurious float32 matmul warnings on arm64
        u = (cr @ t["Ey_re"].astype(dtype) - ci @ t["Ey_im"].astype(dtype)) * t["wy"].astype(dtype)
    if not np.isfinite(u).all():
        raise FloatingPointError("inflow plane is not finite")
    w2 = np.outer(t["wz"], t["wy"]).astype(dtype)
    u[0] -= (u[0].sum(dtype=np.float64) / w2.sum(dtype=np.float64)).astype(dtype) * w2
    return u


def mean_coefficients(t):
    """M_n with plane-mean(u)(time) = Re sum_n M_n exp(i phase_n): lets rows be evaluated without the whole plane."""
    b = (t["B_re"][0].astype(np.float64) + 1j * t["B_im"][0]).sum(0)
    e = ((t["Ey_re"].astype(np.float64) + 1j * t["Ey_im"]) * t["wy"]).sum(1)
    return b * e / (t["wz"].astype(np.float64).sum() * t["wy"].astype(np.float64).sum())


def rows(t, times, layers, ys=None):
    """u'(y, z, time) for a list of times and layers (and optionally a subset of y), float64, shape [nt, 3, nl, ny']."""
    times = np.atleast_1d(np.asarray(times, np.float64))
    ys = np.arange(t["ny"]) if ys is None else np.asarray(ys)
    ph = np.exp(1j * (t["psi"][None, :] - np.outer(times, t["omega"])))    # [nt, n]
    ey = (t["Ey_re"][:, ys].astype(np.float64) + 1j * t["Ey_im"][:, ys])   # [n, ny']
    b = t["B_re"][:, layers].astype(np.float64) + 1j * t["B_im"][:, layers]   # [3, nl, n]
    u = np.einsum("tn,iln,ny->tily", ph, b, ey, optimize=True).real * t["wy"][ys]
    mean = (ph @ mean_coefficients(t)).real                                # [nt]
    u[:, 0] -= mean[:, None, None] * t["wz"][layers][None, :, None] * t["wy"][ys][None, None, :]
    return u


def direct(design, t, time, iy, iz):
    """The formula of the module docstring, written out mode by mode for one cell (without the mean removal):
    an independent check that the stored factors encode it."""
    k, s, q, g = design["k"], design["s"], design["q"], design["g"]
    th = k[:, 1] * design["y"][iy] + k[:, 2] * design["z"][iz] + t["psi"] - t["omega"] * time
    return np.array([t["wy"][iy] * t["wz"][iz] * g[iz, i] * np.sum(q[iz] * s[:, i] * np.cos(th)) for i in range(3)])


def main():
    nb, out = int(sys.argv[1]), Path(sys.argv[2]).resolve()
    sheet = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).resolve().parent.parent / "results" / "tunnel-data.json"
    gain = json.loads(Path(sys.argv[4]).read_text()) if len(sys.argv) > 4 else None
    out.mkdir(parents=True, exist_ok=True)
    tab, d = build(sheet, nb, gain)
    path = out / ("inflow_b%d.bin" % nb)
    write_table(path, tab)
    back = read_table(path)
    assert all(np.array_equal(back[key], tab[key]) for key in ("omega", "psi", "U", "wy", "wz", "B_re", "B_im", "Ey_re", "Ey_im"))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    p, sh = d["prof"], d["sheet"]
    b_m, floor = 32.0, FLOOR_CELLS / nb
    sv = lambda arr, zb: float(arr[int(np.argmin(np.abs(sh["z"] - zb)))]) * b_m   # value at a sheet height, metres
    L_at_sheet = np.maximum(L_FACTOR * d["L_sheet"], floor)
    meta = {"cells_per_b": nb, "ny": tab["ny"], "nz": tab["nz"], "modes": tab["n"], "modes_master": N_MASTER, "seed": SEED,
            "table": path.name, "sha256": sha, "bytes": path.stat().st_size, "gain": gain,
            "units": "lengths in b, speeds in U_H (lattice = value * UREF), time in b/U_H (step * UREF / cells_per_b)",
            "parameters": {"lambda_min_cells": LAMBDA_MIN_CELLS, "k_lo_over_ke_ref": K_LO_OVER_KE, "C_mu": C_MU, "floor_cells": FLOOR_CELLS,
                           "L_factor": L_FACTOR, "slope_half_window": SLOPE_HALF_WINDOW, "U_c": U_C, "taper_cells": TAPER_CELLS,
                           "sweep_rms_U_H": d["us"], "ke_ref_per_b": d["ke_ref"], "k_lo_per_b": d["k_lo"], "k_hi_per_b": d["k_hi"]},
            "sheet": {"z_b": sh["z"].tolist(), "dUdz_UH_per_b": d["dudz_sheet"].tolist(), "L_rule_b": d["L_sheet"].tolist(),
                      "L_rule_m": (d["L_sheet"] * b_m).tolist(), "L_used_b": L_at_sheet.tolist(), "L_used_m": (L_at_sheet * b_m).tolist(),
                      "floor_binds": (L_FACTOR * d["L_sheet"] < floor).tolist()},
            "layers": {"z_b": p["z"].tolist(), "U": p["U"].tolist(), "su": p["su"].tolist(), "sv": p["sv"].tolist(), "sw": p["sw"].tolist(),
                       "L_rule_b": p["L_rule"].tolist(), "L_b": p["L"].tolist(), "L_spec_b": p["L_spec"].tolist(), "floor_binds": p["floor_binds"].tolist(),
                       "ke_per_b": d["ke"].tolist(), "n_eff": d["n_eff"].tolist(), "gain": d["g"].tolist()}}
    if d["u_gain"] is not None:                                            # only with the key "U": without it the sidecar is what it was
        meta["layers"]["U_sheet"], meta["layers"]["U"], meta["layers"]["U_gain"] = p["U"].tolist(), d["u_table"].tolist(), d["u_gain"].tolist()
    (out / ("inflow_b%d.json" % nb)).write_text(json.dumps(meta, indent=1))
    print("table %s: %d modes of %d, %d x %d cells, %.2f MB, sha256 %s" % (path, tab["n"], N_MASTER, tab["ny"], tab["nz"], path.stat().st_size / 1e6, sha[:16]))
    print("length scale by the rule (full-scale analogue, b = 32 m), floor %.1f m = %g cells:" % (floor * b_m, FLOOR_CELLS))
    print("  z/b     z[m]  dU/dz[1/s]  L rule[m]  L used[m]")
    for i, zb in enumerate(sh["z"]):
        print("  %6.4f %6.1f  %9.4f  %8.1f  %8.1f%s" % (zb, zb * b_m, d["dudz_sheet"][i] * 10.0 / b_m, d["L_sheet"][i] * b_m, L_at_sheet[i] * b_m,
                                                        "  floor" if L_FACTOR * d["L_sheet"][i] < floor else ""))
    print("summary: inflow b%d: seed %d, %d modes, table sha256 %s; L by the rule at z/b 0.125 / 2 / 7.5 = %.1f / %.1f / %.1f m, used %.1f / %.1f / %.1f m (floor %.1f m);"
          " U_c = U_H, cut-off wavelength %g cells, taper %d cells" % (nb, SEED, tab["n"], sha[:12], sv(d["L_sheet"], 0.125), sv(d["L_sheet"], 2.0), sv(d["L_sheet"], 7.5),
                                                                      sv(L_at_sheet, 0.125), sv(L_at_sheet, 2.0), sv(L_at_sheet, 7.5), floor * b_m, LAMBDA_MIN_CELLS, TAPER_CELLS))


if __name__ == "__main__":
    main()
