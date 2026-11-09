# Changes against upstream FluidX3D

Base: [FluidX3D](https://github.com/ProjectPhysX/FluidX3D) commit `b85d64ee39a72687d163d42806c9a0ccdaed6d3c` (tag v3.8 plus 10 commits, 2026-10-29). `changes.patch` is the `git diff` from that commit to the files in this repository; `git apply changes.patch` on a checkout of that commit reproduces `src/`.

Four files differ. `makefile`, `make.sh` and every other file under `src/` (including `src/OpenCL/` and `src/X11/`) are byte-identical to upstream. Each altered file starts with a comment line that names the base commit and this study.

| File | Added | Removed | Change |
|---|---|---|---|
| `src/defines.hpp` | 8 | 7 | compile-time switches and frame settings |
| `src/lbm.cpp` | 6 | 2 | colour-scale maximum and Q threshold read from the environment |
| `src/kernel.cpp` | 9 | 4 | single-hue colour scale for vortex surfaces |
| `src/setup.cpp` | 192 | 1363 | the wind-study setup replaces the upstream example setups |
| total | 215 | 1376 | |

## `src/defines.hpp` (+8, -7)

Only `#define` lines changed, plus the header comment; no new macros.

- `BENCHMARK` turned off (upstream ships with it on), so `main_setup()` in `setup.cpp` runs instead of the benchmark.
- `EQUILIBRIUM_BOUNDARIES` turned on: the setup fixes the inflow profile on the inlet, outlet, side and top faces with `TYPE_E` cells.
- `SUBGRID` turned on: the Smagorinsky-Lilly LES model keeps the run stable at the very low lattice viscosity of full-scale air.
- `GRAPHICS` turned on: frames are rendered off-screen and written as PNG; there is no interactive window.
- `GRAPHICS_FRAME_WIDTH` and `GRAPHICS_FRAME_HEIGHT` 1920x1080 to 2560x1440.
- `GRAPHICS_BACKGROUND_COLOR` `0x000000` to `0xEDEDED` (light grey, to match the blog images).

## `src/lbm.cpp` (+6, -2)

- Adds `#include <cstdlib>` and a small helper `fx_env(name, default)` that reads a float from an environment variable.
- In `LBM_Domain::Graphics::device_defines()`, the OpenCL define `def_scale_u` uses `FX_U_MAX` (default `GRAPHICS_U_MAX`, 0.18) and `def_scale_Q_min` uses `FX_Q` (default `GRAPHICS_Q_CRITERION`, 0.0001).

Why: both values are baked into the OpenCL kernels when they are compiled at start-up. Reading them from the environment lets each run set the velocity colour range and the Q-criterion threshold without recompiling the C++ code. `setup.cpp` checks that both parse as numbers and that `FX_U_MAX` is positive before the kernels are built.

## `src/kernel.cpp` (+9, -4)

- Adds the OpenCL function `colorscale_navy(x)`: one hue from light blue (0.80, 0.84, 0.94) to navy `#163795`, with a square-root ramp so slow flow is still visible.
- In the kernels `graphics_q_field` and `graphics_q` (Q-criterion points and marching-cubes isosurfaces), velocity colouring (`field_mode == 0`) uses `colorscale_navy` instead of `colorscale_rainbow`. The rainbow scale still exists and is still used by the other field renderers.

Why: the vortex renderings had to read as one-colour drawings on a light background.

## `src/setup.cpp` (+192, -1363)

The upstream file (the benchmark setup and the commented-out example setups) is replaced by one `main_setup()` for the wind study. Everything it does:

- Reads its parameters from environment variables (list in `README.md`), parses them strictly and checks their ranges; on a violation it prints one line starting with `Error: wind setup:` and exits with code 1. `OUT` gets a trailing slash when it lacks one. Physical units come from one fixed assumption: the reference wind speed at `ZREF_M` is 10 m/s, so one time step is `DX_M * UREF / 10` seconds.
- Loads a 2.5D height map (`HM`, uint16 `HX` x `HY`) and marks every cell with `z < h` as solid: the flat ground layer and the buildings.
- Optionally voxelises a binary STL (`STL`) as solid first, in cell-index coordinates with the origin at the domain corner, and prints its triangle and solid-cell counts.
- Sets a power-law inflow profile `UREF * (max(z - 0.5, 0.5) / ZREF)^ALPHA` (z in cells) in +x in the whole domain and fixes it on the inlet, outlet, side and top faces (`TYPE_E`).
- Runs `T_SI` seconds of wind in `FRAMES` stages. In the second half of the run, every `SAMPLE_SI` seconds (rounded down to whole time steps) it reads back three planes from the GPU: the horizontal slices at `ZPED` and `ZHI` and the vertical section at `y = PY`, and accumulates sums and sums of squares of the speed in double precision.
- At the end of each stage writes raw float32 arrays to `OUT/planes/`: instantaneous speed (`hz_inst-<k>`, `hi_inst-<k>`, `vy_inst-<k>`), `samples-<k>.txt` with the number of samples so far (0 before the averaging window), and, once samples exist, running means (`hz_mean`, `hz_mean_ux`, `hz_mean_uy`, `hi_mean`, `vy_mean`), tagged copies of the means so far (`hz_mean-<k>`, `hi_mean-<k>`) to check that the average has settled, standard deviations of the speed (`hz_std`, `hi_std`) and `meta.txt` with the grid, the unit conversion, the real sampling period, the averaging window and the animation cadence.
- Optionally accumulates an axis-aligned box (`VX0 <= x < VX1`, `VY0 <= y < VY1`, `0 <= z < VZ1`), read back from the device as one contiguous slab per component, and writes its mean velocity (`vol_mean.f32`), last velocity (`vol_inst-<k>.f32`), mean density (`vol_rho_mean.f32`) and `vol_meta.txt`, all in `OUT/planes/`.
- Renders three views per stage with the built-in ray tracer: a close and a wide Q-criterion view (`q3d/`, `q3d_wide/`) and a vertical section at `y = PY` showing the field only (`side/`). Cameras are placed by azimuth, elevation, distance and field of view around the landmark cell `PX, PY`.
- Optionally animates the last `ANIM_T` seconds: every `ANIM_SI` seconds writes the upper slice (`anim/hi-NNNN.f32`) and a frame of the close view (`anim/q-NNNN.png`).

## Fixes after the published results

The first, tagged commit of this repository holds the source exactly as it was compiled for the published images. The current version adds the fixes below, from an independent review (2026-11-09). On the machine that computed the results (an Apple Silicon GPU, whose OpenCL buffers are zero-copy) they change no output except where stated. Diff from the as-used source: 4 files, 65 lines added, 37 removed.

- Box read-back (`setup.cpp`, lambda `read_box`): the box accumulation and `vol_inst` read host copies of `u` and `rho` that, in the as-used version, were only refreshed on the two sampled slices and the `y = PY` row. On a GPU with unified memory the host copy is the device buffer, so the as-used results are correct; on a discrete GPU the box would have held the initial inflow profile. Now the slab `z < VZ1` is read from the device for the three velocity components and for the density before every use; the ranged read is a no-op on zero-copy devices. No change in the outputs here.
- Input validation: every variable is parsed strictly (a value that is not a number, or a negative or fractional value for an integer variable, is an error); `HX`, `HY`, `NZ`, `FRAMES`, `DX_M`, `UREF`, `T_SI` and `FX_U_MAX` must be positive; `ZPED`, `ZHI` below `NZ`; `PX` below `HX`; `PY` below `HY`; when any of `VX0 VX1 VY0 VY1 VZ1` is set, `VX0 < VX1 <= HX`, `VY0 < VY1 <= HY`, `0 < VZ1 <= NZ` must hold (the as-used version silently switched the box off instead). `OUT` gets a trailing slash. A missing height map names the expected byte count. New errors only; valid inputs behave as before.
- Means and standard deviations: the sums of the speed and of its square on both slices are accumulated in double and the standard deviation is computed in double; the as-used version used float32 for both. `hz_mean`, `hi_mean`, `hz_std`, `hi_std` and the tagged means may differ in the last bits. The other means (`hz_mean_ux`, `hz_mean_uy`, `vy_mean`, the box) are unchanged float32 sums.
- `samples-<k>.txt` is written for every stage, with 0 for stages that end before the averaging window; the as-used version wrote it only once samples existed.
- `meta.txt` gets six new lines after the existing ones: `sample_steps`, `sample_period_s`, `avg_start_s`, `avg_end_s`, `anim_start_s`, `anim_period_s` (the last two are 0 when the animation is off). Existing lines are unchanged.
- Wording: the camera target is called the landmark (the as-used version named the building of the first case, the Palace of Culture and Science, in a variable, in comments and in the console line); the comment on the colour scale says navy, not brand navy; each altered file starts with a header line naming the base commit; the comment on the `STL` variable now states the right convention (cell-index coordinates with the origin at the domain corner; the as-used comment said the origin was the domain centre, which the geometry scripts had already found to be wrong and compensated for, so the published geometry was placed correctly).

## Regression test
Regression test (2026-11-09): the as-built city run of the study (750 x 400 x 200 cells of 3 m, 600 s) repeated with the published-results build and with this fixed build on the same Apple M5 Pro: mean-speed planes differ by at most 2.6e-8 lattice units (relative 1.3e-6), the standard-deviation planes by at most 4e-5 (the double accumulators' rounding), the 3D box is byte-identical, and `meta.txt` differs only by the new lines. The box read-back fix is a no-op on a unified-memory GPU, as expected.
