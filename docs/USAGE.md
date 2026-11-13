# Usage

The solver reads everything from environment variables and writes to the folder `OUT`. There is no configuration file and no interactive window.

## Conventions

- Lengths are in metres, angles in degrees, speeds in lattice units unless noted.
- Azimuth is measured from +x towards +y; elevation above the horizon.
- The wind enters from the west, at x = 0, and blows along +x.
- The physical reference speed is fixed in the code at 10 m/s and is used only to convert time and speed: one time step is `DX_M * UREF / 10` seconds, and speeds in lattice units times `10 / UREF` are m/s.

## Errors

The variables added for this study are parsed strictly. An invalid or out-of-range value, a variable that needs another one, or a table of the wrong size prints one line starting with `Error: wind setup:` and stops the solver with exit code 1 before it computes anything.

Four older variables are read without a range check: `SAMPLE_SI`, `NU`, `ZREF_M` and `ALPHA`. Give them sensible positive values.

## Environment variables

### Domain and input

| Variable | Default | Meaning |
|---|---|---|
| `HM` | required | Path of the height map (format below). |
| `HX`, `HY` | 1000, 600 | Height-map size in cells. The domain is `HX` x `HY` x `NZ` cells. |
| `NZ` | 200 | Cells in z. |
| `DX_M` | 3 | Cell size in metres, > 0. |
| `STL` | none | Optional binary STL, voxelised as solid before the height map (format below). |

### Wind

| Variable | Default | Meaning |
|---|---|---|
| `UREF` | 0.05 | Inflow speed at `ZREF_M` in lattice units, > 0. |
| `ZREF_M` | 100 | Reference height of the profile in metres. |
| `ALPHA` | 0.24 | Power-law exponent of the default inflow profile. The profile is `UREF * (max(z - 0.5, 0.5) / ZREF)^ALPHA` with z and `ZREF` in cells: the wall sits at z = 0.5, and the clamp gives the lowest cells the speed of half a cell above the wall. |
| `PROFILE` | none | A tabulated mean profile, used in place of the power law (format below). |
| `INFLOW` | none | A turbulent-inflow table for this grid (format below). The table's mean profile is used on every fixed face and as the initial field (`PROFILE` must be unset; `ALPHA` is ignored). At the plane x = 0 the table's fluctuation is added. Single GPU only. |
| `INFLOW_K` | 10 | Steps between two computed inlet planes; needs `INFLOW`. In between, the inlet velocity is interpolated linearly in time on the GPU while the host computes the next plane. |
| `NU` | 2e-8 | Kinematic viscosity in lattice units. |
| `FX_CS` | upstream value, about 0.1733 | Smagorinsky constant C, 0 < C <= 0.5. Read when the OpenCL kernels are built. |
| `Z0_M` | none | Roughness length of the ground in metres. Required in the `GROUND_WALL_MODEL` build and an error in the other builds. It must be below half a cell. Cd = (0.4 / ln(0.5 DX_M / Z0_M))^2. The first air cell above every column whose height-map value is 1 gets a specular reflection minus the stress `rho Cd \|u\| u`; roofs and walls keep bounce-back. |

### Time, sampling and output

| Variable | Default | Meaning |
|---|---|---|
| `T_SI` | 600 | Simulated wind time in seconds, > 0. The run has `floor(T_SI / dt)` time steps. |
| `SAMPLE_SI` | 1 | Seconds between samples, rounded down to whole time steps, at least one. The real period is written to `meta.txt`. |
| `SPIN_SI` | unset | Seconds of wind before the averaging starts, 0 <= `SPIN_SI` < `T_SI`. Unset: the second half of the run is averaged. Only samples taken after the start of the window enter the means. |
| `FRAMES` | 4 | Number of output stages, > 0. Stage `k` ends at step `steps * k / FRAMES`. At every stage end the run writes the instantaneous planes, `samples-<k>.txt` and `meta.txt` and renders three views. Means, tagged means, standard deviations and second moments are written once at least one sample has been averaged. |
| `OUT` | `export/` next to the binary | Output folder; a missing trailing slash is added. Subfolders `planes/`, `q3d/`, `q3d_wide/`, `side/`, `anim/`. |
| `PX`, `PY` | 633, 300 | Cell of the landmark the cameras look at. The vertical section passes through `y = PY`. `PX < HX`, `PY < HY`. |
| `ZPED` | 1 | z index of the near-ground slice, below `NZ`. |
| `ZHI` | `100 / DX_M`, rounded | z index of the upper slice, below `NZ`. |
| `ZPEDS` | none | `a,b,...`: more near-ground layers, each >= 1 and < `NZ`, none twice. For every layer k the run writes what it writes for `ZPED`, named `hz<k>_mean.f32`, `hz<k>_mean_ux.f32`, `hz<k>_mean_uy.f32`, `hz<k>_std.f32`, `hz<k>_mean-<stage>.f32`; `meta.txt` gets a last line `zpeds_cells`. |
| `XST` | none | `a,b,...`: x-normal planes, cell indices below `HX`. Over the averaging window the run keeps, in double, the sums of u, v, w, uu, vv, ww and uw per cell, and writes them at every stage as `planes/xst-<k>.f64` with `xst_meta-<k>.txt` (per plane seven arrays of `HY * NZ`, y fastest). The planes are gathered on the GPU, so only they are read back. Single GPU only. |

### Cameras

| Variable | Default | Meaning |
|---|---|---|
| `QA`, `QB`, `QD`, `QF`, `QX` | -115, 28, 1500, 40, 60 | Close vortex view: azimuth, elevation, distance, field of view, and target shift along +x from the landmark. The target sits 40 m above ground. Written to `q3d/`. |
| `WA`, `WB`, `WD`, `WF` | -150, 32, 3400, 40 | Wide vortex view, aimed 150 m west of the landmark. Written to `q3d_wide/`. |
| `SD`, `SF` | 1700, 40 | Section view: distance and field of view. It looks north at the plane `y = PY`. Written to `side/`. |
| `RENDER` | 1 | 0 writes no rendered frames (vortex views, section, film images). Every plane is still written. |

### Optional 3D box

| Variable | Default | Meaning |
|---|---|---|
| `VX0`, `VX1`, `VY0`, `VY1`, `VZ1` | none | Cell range `VX0 <= x < VX1`, `VY0 <= y < VY1`, `0 <= z < VZ1`. When any of them is set, all must satisfy `VX0 < VX1 <= HX`, `VY0 < VY1 <= HY`, `0 < VZ1 <= NZ`. The run writes the mean velocity `vol_mean.f32`, the last velocity `vol_inst-<k>.f32`, the mean density `vol_rho_mean.f32` and `vol_meta.txt` (component-major, x fastest, then y, then z). The box is read back from the GPU before every use. |
| `VOL_TAGGED` | 0 | 1 also keeps `vol_mean-<k>.f32` and `vol_meta-<k>.txt` of every stage. |
| `VOL_M2` | 0 | 1 (needs the box) writes `vol_std.f32` at every stage: standard deviations of u, v, w and the covariance of u and w over the samples so far, four arrays of the box, from sums kept in double. |

### Optional film

| Variable | Default | Meaning |
|---|---|---|
| `ANIM_T` | 0 (off) | Length in seconds of the final part of the run to film. |
| `ANIM_SI` | 1 | Seconds between film frames. A film frame is taken only at a sample step, so the real interval is never shorter than `SAMPLE_SI`. `meta.txt` reports the interval rounded to whole solver steps, not the real spacing of the frames. Each frame writes `anim/hi-NNNN.f32` (upper slice speed) and `anim/q-NNNN.png` (close vortex view, unless `RENDER=0`). |
| `ANIM_SECTION` | 0 | 1 (needs `ANIM_T` > 0) also writes u and w on the plane `y = PY` per frame: `anim/sec_u-NNNN.f32` and `anim/sec_w-NNNN.f32` (`HX` x `NZ`, x fastest). |

### Rendering

Both are read when the OpenCL kernels are built at start-up.

| Variable | Default | Meaning |
|---|---|---|
| `FX_U_MAX` | 0.18 | Speed at the top of the colour scale, lattice units, > 0. |
| `FX_Q` | 0.0001 | Q-criterion threshold of the vortex surfaces. |

### Test hooks

Not used for results.

| Variable | Default | Meaning |
|---|---|---|
| `PERIODIC_XY` | 0 | 1: periodic side faces instead of fixed-profile faces. Must be 0 with `INFLOW`. |
| `WIND_DEG` | 0 | The inflow turned about z by this many degrees. Must be 0 with `INFLOW`. |
| `ZPROF_EVERY` | 0 | n > 0: every n steps the sums of density and momentum per layer, written to `planes/zprof.f64`. With it set, the run also takes an averaging sample at every record step in the averaging window. |
| `INFLOW_STEP0` | 0 | The playback step of solver step 0. |
| `INFLOW_DUMP` | none | `a,b,...`: playback steps after which the inlet velocity is read back to `planes/inlet-<step>.f32`. |
| `FLAGS_DUMP` | 0 | 1: the flag field after initialisation, `planes/flags.u8` (one byte per cell; 0x01 solid, 0x02 fixed-profile, 0x40 marked for the wall stress). |

## Input formats

### Height map (`HM`)

A raw array of `HX * HY` unsigned 16-bit little-endian integers, row-major with x fastest (index `y * HX + x`). The value `h` of a column makes the cells `z = 0 ... h - 1` solid: `1` is bare ground (only the ground layer at `z = 0`), `h > 1` is a building that reaches cell `h - 1`.

Every column needs `h >= 1` and `h <= NZ`. This is a requirement on the input: the solver does not check the values or the length of the file. The terrain is flat.

### STL (`STL`)

A binary STL, read at scale 1 without repositioning, in cell-index coordinates: cell `(x, y, z)` is centred at `(x, y, z)` and the origin is the domain corner. The STL is voxelised before the height map; the cells it marks solid stay solid.

The voxeliser casts rays along the axis with the smallest cross-section of the mesh. Along a ray, a cell whose centre lies less than one cell below the last exit face stays air, so a mesh must end at least 0.75 cell beyond the centre of the last cell it is meant to fill.

### Mean profile (`PROFILE`)

A text file of rows `z_m U`: height above the ground in metres, rising, and speed as a fraction of `UREF`. Lines starting with `#` are comments. At least two rows. The profile is linear between rows; outside them the end values hold.

### Turbulent-inflow table (`INFLOW`)

A table of random Fourier modes, made by `validation/scripts/make_inflow_table.py` or `validation/scripts/make_calibrated_inflow.py` for this grid: its plane must be `HY` x `NZ` cells. The format is `write_table()` in `make_inflow_table.py` (little-endian, magic `SCINFLW1`, version 1); the method is in the docstring of that script.

The table is in units of the building width b and of the speed at building height. One solver step is `UREF / nb` of its time unit, where nb is the number of cells across b in the table.

## Outputs

Plane outputs in `OUT/planes/` are raw float32 little-endian arrays of the speed `|u|` in lattice units unless named otherwise.

| File in `OUT/planes/` | Content |
|---|---|
| `hz_*.f32`, `hi_*.f32` | `HX` x `HY` (x fastest) at `ZPED` and `ZHI`. |
| `vy_*.f32` | `HX` x `NZ` at `y = PY`. |
| `hz_mean_ux.f32`, `hz_mean_uy.f32` | Mean components. |
| `hz_std.f32`, `hi_std.f32` | Standard deviation of the speed over the averaging window. |
| `hz<k>_*.f32` | The same for every `ZPEDS` layer k. |
| `xst-<k>.f64`, `xst_meta-<k>.txt` | Sums on the `XST` planes. |
| `vol_*.f32`, `vol_meta*.txt` | The optional 3D box. |
| `samples-<k>.txt` | Number of samples averaged by the end of stage k (0 for stages that end before the averaging window starts). |
| `meta.txt` | The grid, the unit conversion, the slice and landmark indices, the number of samples, the elapsed time, the sampling period, the averaging window and the film cadence. |

Rendered frames go to `q3d/` (close vortex view), `q3d_wide/` (wide vortex view) and `side/` (section), film frames to `anim/`.

## Limits of the setup

- Boundaries: all outer faces (inlet, outlet, both sides and the top) are equilibrium boundaries pinned to the mean inflow profile; there is no outflow condition. Results near the outlet, the sides and the top are confined, so the domain must leave a wide margin around the area of interest.
- Solids: buildings are stair-stepped bounce-back solids on flat ground; terrain is not modelled. The wall-stress condition acts only on the first air cell above open ground.
- The inflow table is made for one grid and one plane size; it has no shear stress and is calibrated at one station.
- Rendered vortex images are not reproducible bit for bit: two runs of one binary on one GPU wrote identical fields and section images, but their vortex images differed in 9 to 140 pixels of 3.7 million. Compare fields, not those images.

The limits of the wind-tunnel check are in [VALIDATION.md](VALIDATION.md#limits).
