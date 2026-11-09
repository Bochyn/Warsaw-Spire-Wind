# Wind around Warsaw Spire, first attempt (altered FluidX3D)

This is an altered copy of [FluidX3D](https://github.com/ProjectPhysX/FluidX3D), the lattice Boltzmann GPU solver by Dr. Moritz Lehmann. An architect used it for an experiment on a personal, non-commercial blog: a qualitative look at the wind around the Warsaw Spire tower and along Łucka street, computed from a height map of the city centre. The images are illustrations, not a wind engineering study. The code is published because the FluidX3D licence requires it (see Licence).

The repository holds two versions. The first, tagged commit is the source exactly as it was compiled for the published results. The current version adds the fixes from an independent review, listed in `CHANGES.md` under "Fixes after the published results"; on the machine that computed the results they change nothing in the outputs apart from rounding in the means and standard deviations, new validation errors and new lines in `meta.txt`.

Regression test (2026-11-09): the as-built city run of the study (750 x 400 x 200 cells of 3 m, 600 s) repeated with the published-results build and with this fixed build on the same Apple M5 Pro: mean-speed planes differ by at most 2.6e-8 lattice units (relative 1.3e-6), the standard-deviation planes by at most 4e-5 (the double accumulators' rounding), the 3D box is byte-identical, and `meta.txt` differs only by the new lines. The box read-back fix is a no-op on a unified-memory GPU, as expected.

## What is altered

Four files differ from upstream commit `b85d64e`; `CHANGES.md` has the details and `changes.patch` the diff.

- `src/setup.cpp`: the upstream example setups are replaced by one setup: a 2.5D city height map (and an optional STL) becomes the solid, a power-law wind profile enters from the west (+x) and is held on the inlet, outlet, side and top faces, the run samples horizontal slices and a vertical section, writes means and standard deviations as raw float32 arrays, optionally a 3D box and an animation, and renders Q-criterion vortex views and a section with the built-in ray tracer. All inputs come from environment variables and are validated.
- `src/defines.hpp`: `BENCHMARK` off; `EQUILIBRIUM_BOUNDARIES`, `SUBGRID` and `GRAPHICS` on; frames 2560x1440 on a light grey background (`#EDEDED`).
- `src/lbm.cpp`: the velocity colour range and the Q-criterion threshold are read at run time from `FX_U_MAX` and `FX_Q`.
- `src/kernel.cpp`: a single-hue colour scale (light blue to navy) replaces the rainbow scale on the Q-criterion vortex renderings.

Everything else, including `makefile`, `make.sh`, `src/OpenCL/` and `src/X11/`, is unchanged.

## Build

The code needs a C++17 compiler and OpenCL. On macOS the OpenCL framework is part of the system. The makefile's `macOS` target also links X11, which this build does not use and which is normally not installed, so the link step fails after all objects have compiled; a failed link removes an existing `bin/FluidX3D`, so keep binaries you still need under another name. Link by hand:

```
make macOS -j$(sysctl -n hw.ncpu)   # compiles src/*.cpp into temp/*.o, then fails at the X11 link step
mkdir -p bin
g++ temp/*.o -o bin/FluidX3D -std=c++17 -pthread -O -framework OpenCL
```

If XQuartz is installed, `make macOS` links on its own. Linux and Android follow the upstream `makefile` and `make.sh` (`make Linux`, `make Linux-X11`, `make Android`). Windows users take the Visual Studio project files from the upstream repository; they are not part of this package.

## Run

The binary reads everything from environment variables and writes to `OUT`. Lengths are in metres, angles in degrees, speeds in lattice units unless noted. Azimuth is measured from +x towards +y; elevation above the horizon. Every variable is parsed strictly; an invalid or out-of-range value prints one line starting with `Error: wind setup:` and exits with code 1.

Domain and input

- `HM`: path of the height map (required).
- `HX`, `HY`: height-map size in cells (default 1000 x 600); the domain is `HX` x `HY` x `NZ` cells.
- `NZ`: cells in z (default 200).
- `DX_M`: cell size in metres (default 3), must be > 0.
- `STL`: optional binary STL voxelised as solid before the height map (default none; format below).

Wind

- `UREF`: inflow speed at `ZREF_M` in lattice units (default 0.05), must be > 0. The physical reference speed is fixed in the code at 10 m/s and is used only to convert time and speed: one time step is `DX_M * UREF / 10` seconds, and speeds in lattice units times `10 / UREF` are m/s.
- `ZREF_M`: reference height of the profile in metres (default 100).
- `ALPHA`: power-law exponent of the inflow profile (default 0.24, terrain category IV of PN-EN 1991-1-4 NA). The profile is `UREF * (max(z - 0.5, 0.5) / ZREF)^ALPHA` with z and `ZREF` in cells: the wall sits at z = 0.5 and the clamp gives the lowest cells the speed of half a cell above the wall instead of zero.
- `NU`: kinematic viscosity in lattice units (default 2e-8; the Smagorinsky model keeps the run stable).

Time, sampling and output

- `T_SI`: simulated wind time in seconds (default 600), must be > 0. The run has `floor(T_SI / dt)` time steps.
- `SAMPLE_SI`: seconds between samples of the planes (default 1). The period is rounded down to whole time steps, at least one: `sample_steps = max(floor(SAMPLE_SI / dt), 1)`; the real period `sample_steps * dt` is written to `meta.txt`. Samples are taken every `sample_steps` steps from the start of each stage and at the stage end (the last interval of a stage can be shorter). Only samples with `t > avg_start_s` are averaged, where `avg_start_s = floor(steps / 2) * dt`: the second half of the run.
- `FRAMES`: number of output stages (default 4), must be > 0. Stage `k` ends at step `steps * k / FRAMES`. At every stage end the run writes the instantaneous planes, `samples-<k>.txt` and `meta.txt` and renders the three views; the means, the tagged means and the standard deviations are written only once at least one sample has been averaged. With `FRAMES=4` the window starts exactly at the end of stage 2, so stages 1 and 2 report 0 samples and the tagged means exist for stages 3 and 4.
- `OUT`: output folder (default `export/` next to the binary; a missing trailing slash is added); subfolders `planes/`, `q3d/`, `q3d_wide/`, `side/`, `anim/`.
- `PX`, `PY`: cell of the landmark the cameras look at; the vertical section passes through `y = PY` (default 633, 300); `PX < HX`, `PY < HY`.
- `ZPED`: z index of the near-ground slice (default 1, the first cell above flat ground); must be below `NZ`.
- `ZHI`: z index of the upper slice (default `100 / DX_M`, rounded); must be below `NZ`.

Cameras

- `QA`, `QB`, `QD`, `QF`, `QX`: close vortex view, azimuth, elevation, distance, field of view and target shift along +x from the landmark (defaults -115, 28, 1500, 40, 60). The target sits 40 m above ground.
- `WA`, `WB`, `WD`, `WF`: wide vortex view (defaults -150, 32, 3400, 40), aimed 150 m west of the landmark.
- `SD`, `SF`: section view, distance and field of view (defaults 1700, 40); it looks north at the plane `y = PY`, aimed 150 m east of the landmark point and 60 m above it (100 m above ground).

Optional 3D box

- `VX0`, `VX1`, `VY0`, `VY1`, `VZ1`: cell range of a box, `VX0 <= x < VX1`, `VY0 <= y < VY1`, `0 <= z < VZ1` (the upper bounds are exclusive). When any of them is set, all must satisfy `VX0 < VX1 <= HX`, `VY0 < VY1 <= HY`, `0 < VZ1 <= NZ`. The run then writes to `OUT/planes/` the mean velocity `vol_mean.f32`, the last velocity `vol_inst-<k>.f32`, the mean density `vol_rho_mean.f32` and `vol_meta.txt` (component-major, x fastest, then y, then z). The box is read back from the device before every use, see Limitations.

Optional animation

- `ANIM_T`: length in seconds of the final part of the run to animate (default 0, off).
- `ANIM_SI`: seconds between animation frames (default 1, rounded down to whole time steps). Each frame writes `anim/hi-NNNN.f32` (upper slice speed) and `anim/q-NNNN.png` (close vortex view).

Rendering (read when the OpenCL kernels are built at start-up)

- `FX_U_MAX`: speed at the top of the colour scale, lattice units (default 0.18), must be > 0.
- `FX_Q`: Q-criterion threshold of the vortex surfaces (default 0.0001).

Height-map format: a raw array of `HX * HY` unsigned 16-bit little-endian integers, row-major with x fastest (index `y * HX + x`). The value `h` of a column makes the cells `z = 0 ... h - 1` solid: `1` is bare ground (only the ground layer at `z = 0`), `h > 1` is a building that reaches cell `h - 1`. Every column needs `h >= 1`; a column with `h = 0` would leave the bottom face without a wall. The terrain is flat; building heights are counted from the ground layer.

STL: a binary STL, read at scale 1 without repositioning, in cell-index coordinates: cell `(x, y, z)` is centred at coordinate `(x, y, z)` and the origin is the domain corner, so the domain spans `-0.5 ... N - 0.5` along each axis. A mesh drawn in domain metres (corner of the domain at 0, cell `i` covering `i * DX_M ... (i + 1) * DX_M`) is divided by `DX_M` and shifted by minus half a cell. The STL is voxelised before the height map; the cells it marks solid stay solid whatever the height map says.

Plane outputs in `OUT/planes/` are raw float32 little-endian arrays of the speed `|u|` in lattice units: `hz_*` and `hi_*` are `Nx` x `Ny` (x fastest) at `ZPED` and `ZHI`; `vy_*` is `Nx` x `Nz` at `y = PY`; `hz_mean_ux`, `hz_mean_uy` are mean components; `hz_std`, `hi_std` the standard deviation of the speed over the averaging window. `meta.txt` holds the grid (`nx ny nz dx_m`), the conversion (`u_ref_lattice u_ref_si z_ref_m`), the slice and landmark indices, the number of samples, the elapsed time `t_si`, the real sampling period (`sample_steps`, `sample_period_s`), the averaging window (`avg_start_s`, `avg_end_s`) and the animation cadence (`anim_start_s`, `anim_period_s`, both 0 when off).

Example (the 3 m city run, two output stages, animation over the last 120 s):

```
cd bin
HM=../heightmap.u16 HX=1000 HY=600 NZ=200 DX_M=3 PX=633 PY=300 T_SI=720 FRAMES=2 SAMPLE_SI=1 \
ANIM_T=120 ANIM_SI=0.5 FX_Q=0.0005 FX_U_MAX=0.15 OUT=../out/ ./FluidX3D
```

## Limitations

- Boundaries: all outer faces (inlet, outlet, both sides and the top) are equilibrium boundaries pinned to the inflow profile; there is no outflow condition. Results near the outlet, the sides and the top are confined and untrustworthy, so the domain must leave a wide margin around the area of interest. Our runs left 1.1 km behind the tower at 3 m cells and 550 m at 1 m cells.
- Inflow: a steady power-law profile without synthetic turbulence and without surface roughness.
- Solids: buildings are stair-stepped bounce-back solids on flat ground; terrain is not modelled.
- Averaging: the means cover the second half of the run. The standard deviations measure the resolved unsteadiness of a smooth inflow (vortex shedding, wakes), not atmospheric gusts.
- No wall model: near-ground speeds and surface forces are uncertain.
- Box read-back: the current version reads the 3D box from the device before using it. Versions before this fix (including the one that produced the published results) read host memory that only the sampled slices refreshed; they were correct only on GPUs with unified memory, where the host copy is the device buffer, and would have written the initial inflow profile into the box on a discrete GPU.

## References

The FluidX3D licence asks for these publications to be cited when results are published; the blog post cites them. The list is copied from the upstream README.

- Lehmann, M.: [Computational study of microplastic transport at the water-air interface with a memory-optimized lattice Boltzmann method](https://doi.org/10.15495/EPub_UBT_00006977). PhD thesis, (2023)
- Lehmann, M.: [Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs](https://doi.org/10.3390/computation10060092). Computation, 10, 92, (2022)
- Lehmann, M., Krause, M., Amati, G., Sega, M., Harting, J. and Gekle, S.: [Accuracy and performance of the lattice Boltzmann method with 64-bit, 32-bit, and customized 16-bit number formats](https://www.researchgate.net/publication/362275548_Accuracy_and_performance_of_the_lattice_Boltzmann_method_with_64-bit_32-bit_and_customized_16-bit_number_formats). Phys. Rev. E 106, 015308, (2022)
- Lehmann, M.: [Combined scientific CFD simulation and interactive raytracing with OpenCL](https://www.researchgate.net/publication/360501260_Combined_scientific_CFD_simulation_and_interactive_raytracing_with_OpenCL). IWOCL'22: International Workshop on OpenCL, 3, 1-2, (2022)
- Lehmann, M., Oehlschlägel, L.M., Häusl, F., Held, A. and Gekle, S.: [Ejection of marine microplastics by raindrops: a computational and experimental study](https://doi.org/10.1186/s43591-021-00018-8). Micropl.&Nanopl. 1, 18, (2021)
- Lehmann, M.: [High Performance Free Surface LBM on GPUs](https://doi.org/10.15495/EPub_UBT_00005400). Master's thesis, (2019)
- Lehmann, M. and Gekle, S.: [Analytic Solution to the Piecewise Linear Interface Construction Problem and Its Application in Curvature Calculation for Volume-of-Fluid Simulation Codes](https://doi.org/10.3390/computation10020021). Computation, 10, 21, (2022)

## Licence

This is an altered version of FluidX3D by Dr. Moritz Lehmann, https://github.com/ProjectPhysX/FluidX3D, based on commit `b85d64ee39a72687d163d42806c9a0ccdaed6d3c` (v3.8 plus 10 commits). It is not the original software and is marked as altered here, in `CHANGES.md` and in the header line of every altered file (clause 1 of the licence). This repository exists to satisfy clause 5: results generated with an altered version were published, so the altered source is published too.

The use is non-commercial (clause 2): a personal blog, no product, no service, no fee.

`LICENSE.md` is the upstream licence file, unchanged, and applies to this repository (clause 7).

## Not included

- The city height maps, built from Polish GUGiK elevation data and building footprints.
- The Rhino and plotting scripts that produced the published images and the STL of the tower.
- The results (raw arrays, frames, animations).
