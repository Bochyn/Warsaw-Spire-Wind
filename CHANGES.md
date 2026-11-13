# Changes against upstream FluidX3D

Base: [FluidX3D](https://github.com/ProjectPhysX/FluidX3D) commit `b85d64ee39a72687d163d42806c9a0ccdaed6d3c` (tag v3.8 plus 10 commits, 2026-10-29). `changes.patch` is the diff from that commit to the files of this repository: `git apply changes.patch` on a checkout of that commit reproduces `src/` and adds `build.sh`.

Five files under `src/` differ from upstream and one file is new. Every other file under `src/` (including `src/OpenCL/` and `src/X11/`) is byte-identical to upstream. Each altered file starts with a comment line that names the base commit and the study. The upstream `makefile` and `make.sh` are not part of this repository and were not used: `build.sh` replaces them.

| File | Added | Removed | Change |
|---|---|---|---|
| `src/defines.hpp` | 15 | 7 | compile-time switches, frame settings, two build flags |
| `src/lbm.hpp` | 14 | 0 | declarations for the inflow playback, the plane recording and the wall-stress ground |
| `src/lbm.cpp` | 37 | 2 | run-time colour range, Q threshold and Smagorinsky constant; inflow playback and plane recording |
| `src/kernel.cpp` | 51 | 5 | single-hue colour scale; wall-stress ground; two new kernels |
| `src/setup.cpp` | 504 | 1346 | the wind-study setup replaces the upstream example setups |
| `build.sh` | 33 | new | build script |
| total | 654 | 1360 | |

## `src/defines.hpp` (+15, -7)

- `BENCHMARK` turned off (upstream ships with it on), so `main_setup()` in `setup.cpp` runs instead of the benchmark.
- `EQUILIBRIUM_BOUNDARIES` turned on: the setup fixes the inflow profile on the inlet, outlet, side and top faces with `TYPE_E` cells.
- `SUBGRID` turned on: the Smagorinsky-Lilly LES model keeps the run stable at the very low lattice viscosity of full-scale air.
- `GRAPHICS` turned on: frames are rendered off-screen and written as PNG; there is no interactive window.
- `GRAPHICS_FRAME_WIDTH` and `GRAPHICS_FRAME_HEIGHT` from 1920x1080 to 2560x1440.
- `GRAPHICS_BACKGROUND_COLOR` from `0x000000` to `0xEDEDED` (light grey).
- `#define FP16S` is wrapped in `#ifndef FP32_STORAGE`: the compiler flag `-DFP32_STORAGE` keeps the distribution functions in 32-bit floats (control runs and the unit test of the wall-stress ground). Without the flag the storage is FP16S, as before.
- A commented-out switch `GROUND_WALL_MODEL`, off by default; it is set with the compiler flag `-DGROUND_WALL_MODEL`.
- A guard: `GROUND_WALL_MODEL` stops the build unless the velocity set is D3Q19.

## `src/lbm.hpp` (+14, -0)

This file had no header comment in the source as it was built and run; the header line was added in this repository to mark it as altered.

- `extern float fx_ground_cd` (only with `GROUND_WALL_MODEL`): the drag coefficient of the wall-stress ground, set by the setup.
- In `LBM_Domain`: the kernels `kernel_set_inlet` and `kernel_gather_planes`; the buffers `inlet` (two inlet planes, three components each), `xst_buf` and `xst_x` (velocity on listed x-normal planes and their x indices); the counters `inlet_k` and `inlet_step0`; the functions `enable_inlet()`, `enable_xst()` and `gather_planes()`.

## `src/lbm.cpp` (+37, -2)

- Adds `#include <cstdlib>` and a helper `fx_env(name, default)` that reads a float from an environment variable.
- In `LBM_Domain::Graphics::device_defines()`, the OpenCL define `def_scale_u` uses `FX_U_MAX` (default `GRAPHICS_U_MAX`, 0.18) and `def_scale_Q_min` uses `FX_Q` (default `GRAPHICS_Q_CRITERION`, 0.0001). Both values are baked into the OpenCL kernels when they are compiled at start-up; reading them from the environment lets each run set them without recompiling.
- Defines `fx_ground_cd` (only with `GROUND_WALL_MODEL`).
- `enable_inlet(k, step0)` allocates the two-plane inlet buffer (never zero-copy, so a write to the device is a real transfer on every GPU) and creates the kernel `set_inlet`. When the playback is on, `enqueue_stream_collide()` runs `set_inlet` before every time step, with the share of the way between the two stored planes.
- `enable_xst(xs)` and `gather_planes()`: the velocity on the listed x-normal planes is copied on the device into a compact buffer, and only that buffer is read back to the host.
- In the OpenCL defines: `def_smag`, the Smagorinsky coefficient 18·sqrt(2)·C², from the environment variable `FX_CS` (C); without `FX_CS` it is the upstream literal 0.76421222 (C = 0.1733). With `GROUND_WALL_MODEL`: the defines `GROUND_WALL_MODEL` and `def_wm_cd`.

## `src/kernel.cpp` (+51, -5)

- Adds the OpenCL function `colorscale_navy(x)`: one hue from light blue (0.80, 0.84, 0.94) to navy `#163795`, with a square-root ramp.
- In the kernels `graphics_q_field` and `graphics_q` (Q-criterion points and isosurfaces), velocity colouring uses `colorscale_navy` instead of `colorscale_rainbow`. The other renderers keep the rainbow scale.
- The Smagorinsky literal `0.76421222f` in the collision is replaced by `def_smag` (same value unless `FX_CS` is set).
- `apply_ground_stress()` (only with `GROUND_WALL_MODEL`, D3Q19): for a cell of the first fluid layer above open ground (flag `TYPE_X`, set by the setup), the four bounced-back diagonal populations are swapped into a local specular reflection, which returns the tangential momentum, and then the log-law wall stress `rho·Cd·|u|·u` is taken out of the horizontal momentum. Mass and vertical momentum of the cell are unchanged. It is called in `stream_collide` and `update_fields`. Roofs and walls keep bounce-back.
- New kernel `set_inlet`: writes the velocity of the fixed cells of the plane x = 0, interpolated linearly in time between the two stored inlet planes. Solid cells of that plane are left alone.
- New kernel `gather_planes`: copies the velocity on the listed x-normal planes into a compact buffer.

## `src/setup.cpp` (+504, -1346)

The upstream file (the benchmark setup and the commented-out example setups) is replaced by one `main_setup()` for the wind study. Its core:

- Reads its parameters from environment variables (list in [docs/USAGE.md](docs/USAGE.md)), parses them strictly and checks their ranges; on a violation it prints one line starting with `Error: wind setup:` and exits with code 1. Physical units come from one fixed assumption: the reference wind speed at `ZREF_M` is 10 m/s, so one time step is `DX_M * UREF / 10` seconds.
- Loads a 2.5D height map (`HM`) and marks every cell with `z < h` as solid: the flat ground layer and the buildings. Optionally voxelises a binary STL (`STL`) as solid first.
- Sets a power-law inflow profile in +x in the whole domain and fixes it on the inlet, outlet, side and top faces (`TYPE_E`).
- Runs `T_SI` seconds of wind in `FRAMES` stages, samples two horizontal slices and the vertical section at `y = PY`, accumulates means and standard deviations of the speed in double precision, and writes them as raw float32 arrays with `meta.txt`.
- Optionally accumulates a 3D box (`VX0` to `VZ1`), read back from the device, and writes its mean velocity, last velocity, mean density and `vol_meta.txt`.
- Renders a close and a wide Q-criterion view and a vertical section per stage with the built-in ray tracer, and optionally a film of the last `ANIM_T` seconds.

Options, all off unless their variable is set (meanings in [docs/USAGE.md](docs/USAGE.md)):

- Turbulent inflow (`INFLOW`, `INFLOW_K`): reads a table of random Fourier modes (format: the writer `validation/scripts/make_inflow_table.py`; the reader here follows `validation/scripts/inflow_reference.cpp`). The table's mean profile is used on every fixed face and as the initial field. At the plane x = 0 the fluctuation of the table is added: every `INFLOW_K` steps the host computes the next inlet plane in a second thread (rows on up to 8 threads), uploads it, and the kernel `set_inlet` interpolates between two planes in time. The mean of the fluctuation over the plane is removed, so the net inflow does not fluctuate. The run is cut into parts that end at every plane update, so samples are taken at the same steps as without `INFLOW`. Test hooks `INFLOW_STEP0` and `INFLOW_DUMP`.
- Tabulated mean profile (`PROFILE`): rows of height in metres and speed as a fraction of `UREF`, linear between rows, end values outside; used in place of the power law.
- Wall-stress ground (build flag `GROUND_WALL_MODEL`, `Z0_M`): Cd = (0.4 / ln(0.5 DX_M / Z0_M))²; the first air cell above every column whose height-map value is 1 (open ground) is flagged `TYPE_X` for the kernel. The setup counts the flagged cells and stops if one of them has a non-solid cell below or beside it in the ground layer. Without the build flag, setting `Z0_M` is an error.
- Second moments on x-normal planes (`XST`): over the averaging window the sums of u, v, w, uu, vv, ww, uw per cell are kept in double and written at every stage (`xst-<k>.f64`, `xst_meta-<k>.txt`).
- In the 3D box: `VOL_TAGGED` keeps the box mean of every stage; `VOL_M2` writes standard deviations of u, v, w and the u-w covariance (`vol_std.f32`) from sums kept in double.
- `SPIN_SI`: seconds of wind before the averaging starts (default: the second half of the run is averaged, as before).
- `ZPEDS`: the near-ground outputs of `ZPED` written for more layers (`hz<k>_mean.f32` and the rest); `meta.txt` gets a line `zpeds_cells`.
- `RENDER=0` switches off all rendered frames; `ANIM_SECTION=1` writes u and w on the plane `y = PY` for every film frame.
- Test hooks `PERIODIC_XY`, `WIND_DEG`, `ZPROF_EVERY`, `FLAGS_DUMP`.
- Printed at the end with `INFLOW`: the number of inlet planes uploaded and the time the solver waited for them.

The output-only changes (`SPIN_SI`, `ZPEDS`, `RENDER`, `ANIM_SECTION`, `FLAGS_DUMP`, the device-side gathering of the `XST` planes) were checked before they were used: the binary of the source before them and the binary after them ran the same 60 s turbulent case with a block on the same GPU, and every numerical output was identical byte for byte.

## `build.sh` (new)

Compiles `src/*.cpp` into `temp/<variant>/` and links `bin/FluidX3D-wind[-<variant>]`; extra arguments go to the compiler (`./build.sh wall-model -DGROUND_WALL_MODEL`). It never starts the solver. On macOS it links the OpenCL framework; on Linux `-lOpenCL`. Upstream's `make.sh` starts the binary after building it and overwrites `bin/FluidX3D`; this script does neither.

## Record of the build that produced the results

- Every published result of this version was computed with the wall-model build. It was made with the command `./build.sh wm -DGROUND_WALL_MODEL`, which then named the binary `bin/FluidX3D_v2_wm`, using g++ 13.3 on Ubuntu 24.04. That binary has sha256 `fd50b633b7c44e11aec764bff36455042fac5464aec5f71297ebb8a75e73bcd9`. Where the build was repeated on another machine, it gave the same binary byte for byte.
- For this repository only the binary names were changed. The same build is now `./build.sh wall-model -DGROUND_WALL_MODEL`, which writes `bin/FluidX3D-wind-wall-model`. The compiler flags are the same; only the names of the binary and of its object folder differ.
- For the wind-tunnel runs the run records name the hash of the source that was built. For the city runs they name the uploaded binary by its checksum, which is the binary of the wind-tunnel runs; no per-run source hash was recorded there.

## Differences between `src/` and the source that was built

`src/` differs from the source that was built only in comments:

- `src/lbm.hpp` was altered but had no header line saying so; the line was added for this repository.
- Sixteen comment lines that pointed to internal working files, or used internal names of the validation steps, were reworded so that they point to files of this repository: 5 in `src/setup.cpp`, 4 in `src/lbm.hpp`, 3 in `src/kernel.cpp`, 2 in `src/lbm.cpp` and 2 in `src/defines.hpp`. No line was added or removed by this.

Check: with every comment removed, the five altered files are identical to the source that was built (`src/lbm.hpp` apart from the empty first line left by its header). Both sources, compiled on one machine (Apple clang 21) with the flags of `build.sh` in all three variants, give byte-identical object files.

## History

The setup was made in two steps.

- The earlier version (7 to 9 November 2026) had a smooth power-law inflow, the bounce-back ground of upstream and the built-in Smagorinsky constant. It contained the compile-time switches and frame settings of `src/defines.hpp` (without `FP32_STORAGE` and `GROUND_WALL_MODEL`), the colour range and Q threshold read at run time in `src/lbm.cpp`, the single-hue colour scale in `src/kernel.cpp`, and the core of `src/setup.cpp`. It scored 0.57 against a pass mark of 0.66 in the wind-tunnel check ([docs/VALIDATION.md](docs/VALIDATION.md)). The blog post shows its wind-tunnel score as the first attempt, so its source is published too: the tag `v1-first-attempt` of this repository (the source as compiled, then the fixes after an independent review and their regression test).
- This version (10 to 12 November 2026) added everything else: the turbulent inflow, the wall-stress ground, the Smagorinsky constant read at run time, the 32-bit storage flag, the options of `src/setup.cpp`, all of the changes in `src/lbm.hpp`, and `build.sh`.

With none of the options of this version set and without the `GROUND_WALL_MODEL` build flag, the solver is meant to compute what the earlier version computed. This was checked byte for byte on a small case twice: against a build of the earlier version's source, and, on the GPU type of the published runs, against a build of this source from before the turbulent inflow was added.
