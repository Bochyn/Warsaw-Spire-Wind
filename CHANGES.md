# Changes against upstream FluidX3D

Base: [FluidX3D](https://github.com/ProjectPhysX/FluidX3D) commit `b85d64ee39a72687d163d42806c9a0ccdaed6d3c` (tag v3.8 plus 10 commits, 2026-10-29). `changes.patch` is the `git diff` from that commit to the files in this repository; `git apply changes.patch` on a checkout of that commit reproduces `src/`.

Four files differ. `makefile`, `make.sh` and every other file under `src/` (including `src/OpenCL/` and `src/X11/`) are byte-identical to upstream.

| File | Added | Removed | Change |
|---|---|---|---|
| `src/defines.hpp` | 7 | 7 | compile-time switches and frame settings |
| `src/lbm.cpp` | 5 | 2 | colour-scale maximum and Q threshold read from the environment |
| `src/kernel.cpp` | 8 | 4 | single-hue colour scale for vortex surfaces |
| `src/setup.cpp` | 168 | 1364 | the wind-study setup replaces the upstream example setups |
| total | 188 | 1377 | |

## `src/defines.hpp` (+7, -7)

Only `#define` lines changed; no new macros.

- `BENCHMARK` turned off (upstream ships with it on), so `main_setup()` in `setup.cpp` runs instead of the benchmark.
- `EQUILIBRIUM_BOUNDARIES` turned on: the setup fixes the inflow profile on the inlet, outlet, side and top faces with `TYPE_E` cells.
- `SUBGRID` turned on: the Smagorinsky-Lilly LES model keeps the run stable at the very low lattice viscosity of full-scale air.
- `GRAPHICS` turned on: frames are rendered off-screen and written as PNG; there is no interactive window.
- `GRAPHICS_FRAME_WIDTH` and `GRAPHICS_FRAME_HEIGHT` 1920x1080 to 2560x1440.
- `GRAPHICS_BACKGROUND_COLOR` `0x000000` to `0xEDEDED` (light grey, to match the blog images).

## `src/lbm.cpp` (+5, -2)

- Adds `#include <cstdlib>` and a small helper `fx_env(name, default)` that reads a float from an environment variable.
- In `LBM_Domain::Graphics::device_defines()`, the OpenCL define `def_scale_u` uses `FX_U_MAX` (default `GRAPHICS_U_MAX`, 0.18) and `def_scale_Q_min` uses `FX_Q` (default `GRAPHICS_Q_CRITERION`, 0.0001).

Why: both values are baked into the OpenCL kernels when they are compiled at start-up. Reading them from the environment lets each run set the velocity colour range and the Q-criterion threshold without recompiling the C++ code.

## `src/kernel.cpp` (+8, -4)

- Adds the OpenCL function `colorscale_navy(x)`: one hue from light blue (0.80, 0.84, 0.94) to navy `#163795`, with a square-root ramp so slow flow is still visible.
- In the kernels `graphics_q_field` and `graphics_q` (Q-criterion points and marching-cubes isosurfaces), velocity colouring (`field_mode == 0`) uses `colorscale_navy` instead of `colorscale_rainbow`. The rainbow scale still exists and is still used by the other field renderers.

Why: the vortex renderings had to read as one-colour drawings on a light background.

## `src/setup.cpp` (+168, -1364)

The upstream file (the benchmark setup and the commented-out example setups) is replaced by one `main_setup()` for the wind study. Everything it does:

- Reads its parameters from environment variables (list in `README.md`), without validating them. Physical units come from one fixed assumption: the reference wind speed at `ZREF_M` is 10 m/s, so one time step is `DX_M * UREF / 10` seconds.
- Loads a 2.5D height map (`HM`, uint16 `HX` x `HY`) and marks every cell with `z < h` as solid: the flat ground layer and the buildings.
- Optionally voxelises a binary STL (`STL`) as solid first and prints its triangle and solid-cell counts. The voxeliser works in cell-index coordinates with the origin at the domain corner; the comment next to the variable in this version says the origin is the domain centre, which is wrong (the geometry for the published results was built with the corner convention, see `README.md`).
- Sets a power-law inflow profile `UREF * (max(z - 0.5, 0.5) / ZREF)^ALPHA` (z in cells) in +x in the whole domain and fixes it on the inlet, outlet, side and top faces (`TYPE_E`).
- Runs `T_SI` seconds of wind in `FRAMES` stages. In the second half of the run, every `SAMPLE_SI` seconds (rounded down to whole time steps) it reads back three planes from the GPU: the horizontal slices at `ZPED` and `ZHI` and the vertical section at `y = PY`, and accumulates sums and sums of squares of the speed in float32.
- At the end of each stage writes raw float32 arrays to `OUT/planes/`: instantaneous speed (`hz_inst-<k>`, `hi_inst-<k>`, `vy_inst-<k>`) and `meta.txt` with the grid and the unit conversion; once samples exist, also running means (`hz_mean`, `hz_mean_ux`, `hz_mean_uy`, `hi_mean`, `vy_mean`), tagged copies of the means so far (`hz_mean-<k>`, `hi_mean-<k>`, with `samples-<k>.txt`) to check that the average has settled, and standard deviations of the speed (`hz_std`, `hi_std`).
- Optionally accumulates an axis-aligned box (`VX0 <= x < VX1`, `VY0 <= y < VY1`, `0 <= z < VZ1`) from the host copies of the velocity and density (correct only on GPUs with unified memory, see Limitations in `README.md`) and writes its mean velocity (`vol_mean.f32`), last velocity (`vol_inst-<k>.f32`), mean density (`vol_rho_mean.f32`) and `vol_meta.txt`, all in `OUT/planes/`.
- Renders three views per stage with the built-in ray tracer: a close and a wide Q-criterion view (`q3d/`, `q3d_wide/`) and a vertical section at `y = PY` showing the field only (`side/`). Cameras are placed by azimuth, elevation, distance and field of view around the landmark cell `PX, PY`.
- Optionally animates the last `ANIM_T` seconds: every `ANIM_SI` seconds writes the upper slice (`anim/hi-NNNN.f32`) and a frame of the close view (`anim/q-NNNN.png`).

Comments and console messages say "PKiN preview": the Palace of Culture and Science was the first case run with this setup; the Warsaw Spire study reused it with another height map and an STL of the tower.

## Fixes after the published results

None yet. This is the source exactly as compiled for the published images.
