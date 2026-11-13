# Validation against a wind-tunnel benchmark

## Summary

The altered solver was checked against AIJ benchmark case A (Tominaga et al. 2008; tunnel data Meng and Hibi 1998): a single block of proportions 1:1:2 in a turbulent boundary layer, wind normal to a face. The tunnel measured the mean velocity at 181 locations in the vertical centre plane and in two horizontal planes, 0.125 and 1.25 building widths above the ground. The simulation is a full-scale analogue with b = 32 m.

Read as a whole: the single block passes from 24 cells across it; at 16 cells four of five rows are met. A second case, AIJ case D (a tower among city blocks), was scored once and did not pass. The full list of limits is at the end of this page.

## Terms

| Term | Meaning |
|---|---|
| b | The width of the block. The block is b x b x 2b; at full scale b = 32 m. Heights and distances are given in b, for example z/b. |
| Approach speed, U_H | The mean wind speed at building height (z = 2b) in the approach flow. Speeds are divided by it. |
| Cells across b | The grid: 16, 24 or 32 cells across the block width, that is cells of 2, 1.333 and 1 m. |
| Location | A point where the tunnel measured the mean velocity. There are 181. |
| Hit | The simulated value is within 25% of the tunnel value or within 0.05 of the approach speed (the COST 732 criterion). |
| Hit rate | The share of locations that are hits. The interval in brackets is the Wilson 95% interval. |
| Pedestrian level | z/b = 0.125, 4 m at full scale. |
| Regions | Upstream (x/b < -0.5), sides, roof, wake below roof height, shear layer above the wake, and pedestrian level, which overlaps the others. |
| Reattachment length | Where the mean streamwise wind on the centre line behind the block turns positive again, measured from the leeward face. |
| Empty domain | The same domain without the block, used to calibrate and judge the approach flow. |
| Inflow sheet | The sheet of the tunnel workbook with the approach-flow profile: mean speed, three standard deviations and k at 24 heights. |
| Pass mark | The threshold a row must meet. Rows A1 to A3, B1 to B6 and C1 below. |
| Declared | Written down before the run it applies to. |
| Earlier version | The first version of this setup: smooth inflow, bounce-back ground, built-in Smagorinsky constant. See [CHANGES.md](../CHANGES.md#history). |

## Pass marks

| Row | What is judged | Pass mark |
|---|---|---|
| A1 | Mean speed of the empty-domain approach flow at the block position, against the inflow sheet | within 5% from z/b = 0.125 up, within 10% at z/b = 0.0625 |
| A2 | Turbulence intensity at building height at the block position | 18% ± 3 points |
| A3 | Change of the mean speed along the last two building heights of the fetch (x/b = -4.5 to -0.5), at every height of the sheet | below 5% |
| B1 | Hit rate of the horizontal speed at all 181 locations | at least 0.66 |
| B2 | Hit rate at pedestrian level | at least 0.66 |
| B3 | Hit rate behind the block below roof height, on the signed streamwise component | at least 0.60 |
| B4 | Hit rate of the weakest region | no region below 0.50 |
| B5 | Reattachment length behind the block at pedestrian level | 0.99 to 1.65 b |
| B6 | Settling of the mean between the last two quarters of the averaging window | mean difference at most 0.01 and largest at most 0.05 of the approach speed (as written; see the results) |
| C1 | Change of the reattachment length on two successive grids | below 10% |

The pass marks of rows B1 to B5 and C1 were fixed before any run of this version. Rows A1 to A3 judge the approach flow in the empty domain; rows B1 to B6 judge a run with the block.

## Data

The Architectural Institute of Japan publishes the benchmark workbooks of its CFD guideline at https://www.aij.or.jp/jpn/publish/cfdguide/index_e.htm. Case A is the single 1:1:2 block; our copy of its workbook, `CaseA_1_1_2.xls`, is 492 544 bytes, sha256 `c2fab8c1494e777381c3819cb8874a10d3f1c5595f701d3f83119d09bad74cff`. Case D, the second check, is the workbook of the high-rise among city blocks on the same page.

The workbooks are not part of this repository.

## What is in `validation/`

| File | What it does |
|---|---|
| `scripts/read_tunnel_data.py` | Reads the case A workbook into JSON: 186 measured rows at 181 locations, and the inflow sheet. |
| `scripts/make_block_geometry.py` | Writes the height map and grid of the block for a given number of cells across b. |
| `scripts/fit_ground_roughness.py` | Fits the roughness length of the ground to the lowest three heights of the inflow sheet. |
| `scripts/make_inflow_table.py` | Makes a turbulent-inflow table of random Fourier modes; holds the file format and a reference evaluator in Python. |
| `scripts/make_calibrated_inflow.py` | Makes the same table with a calibrated length scale. It made the tables of the check. |
| `scripts/inflow_reference.cpp` | The standalone C++ reference of the inflow playback that `src/setup.cpp` follows. |
| `scripts/score_block_run.py` | Scores a run with the block: rows B1 to B6. |
| `scripts/judge_approach_flow.py` | Judges an empty-domain run: rows A1 to A3. It can also write the next calibration file. |
| `scripts/report_approach_profiles.py` | Prints and writes the approach-flow profiles of any run with x-normal planes. |
| `run-settings/*.env` | The five run settings of the check. |
| `inflow-calibration.json` | The calibration of the turbulent inflow: multipliers by height of the inflow sheet's three standard deviations and of its mean speed. |

The docstring of each script describes its inputs and rules.

## Repeating the check

You need Python 3 with numpy, and xlrd for the first step; the commands use [uv](https://docs.astral.sh/uv/). The runs need a GPU with OpenCL. Run every command from `validation/`, with the downloaded workbook saved there as `CaseA_1_1_2.xls`.

### Step 1: prepare the inputs (CPU, a few minutes)

```
mkdir -p results inputs runs
# 1. the workbook -> 186 measured rows at 181 locations, and the inflow sheet
uv run --with xlrd python -I scripts/read_tunnel_data.py CaseA_1_1_2.xls results/tunnel-data.json
# 2. height maps and grids at 16, 24 and 32 cells across the block (b = 32 m at full scale)
for n in 16 24 32; do
  uv run --with numpy python -I scripts/make_block_geometry.py $n inputs/block-$n-cells
done
# 3. the empty domains at 16 and 24 cells: the same height map with every value 1, and the same meta.json
for n in 16 24; do
  mkdir -p inputs/empty-domain-$n-cells
  cp inputs/block-$n-cells/meta.json inputs/empty-domain-$n-cells/
  uv run --with numpy python -I -c "import numpy as np, sys; np.ones_like(np.fromfile(sys.argv[1], np.uint16)).tofile(sys.argv[2])" \
    inputs/block-$n-cells/height.u16 inputs/empty-domain-$n-cells/height.u16
done
# 4. the ground roughness from the lowest three heights of the inflow sheet (the run settings carry it as Z0_M)
uv run python -I scripts/fit_ground_roughness.py results/tunnel-data.json results/ground-roughness.json
# 5. the turbulent-inflow tables: calibration file inflow-calibration.json, length-scale factor 2.6
for n in 16 24 32; do
  uv run --with numpy python -I scripts/make_calibrated_inflow.py $n inputs/inflow results/tunnel-data.json inflow-calibration.json 2.6
done
```

Expected results:

| Output | Expected |
|---|---|
| `results/ground-roughness.json` | z0 = 2.947e-4 m |
| `inputs/inflow/inflow_b16.bin` | sha256 `d290a31b8d444328466315a10d8cc9fc43c1625d0fccd912f16ed28d713d0d91` |
| `inputs/inflow/inflow_b24.bin` | sha256 `8d3ac250290d14fd4d2d1301cd34712766102323cead3117f1935954c934398b` |
| `inputs/inflow/inflow_b32.bin` | sha256 `76e5b0d7c2e79f5e9fd7f85f625f37fff48857cb34e6d81b46af20a0243c26ce` |

These are the tables of the scored runs. Repeating step 5 gave the same hashes, but another maths library may change a last bit of a stored float, so compare them.

### Step 2: build and run (GPU)

```
(cd .. && ./build.sh wall-model -DGROUND_WALL_MODEL)
for s in block-16-cells block-24-cells block-32-cells empty-domain-16-cells empty-domain-24-cells; do
  (set -a; . run-settings/$s.env; set +a; OUT=runs/$s/ ../bin/FluidX3D-wind-wall-model)
done
```

| Run setting | Geometry | Inflow table | Judged by |
|---|---|---|---|
| `block-16-cells.env` | `inputs/block-16-cells` | `inflow_b16.bin` | `score_block_run.py` (rows B1 to B6) |
| `block-24-cells.env` | `inputs/block-24-cells` | `inflow_b24.bin` | `score_block_run.py` |
| `block-32-cells.env` | `inputs/block-32-cells` | `inflow_b32.bin` | `score_block_run.py` |
| `empty-domain-16-cells.env` | `inputs/empty-domain-16-cells` | `inflow_b16.bin` | `judge_approach_flow.py` (rows A1 to A3) |
| `empty-domain-24-cells.env` | `inputs/empty-domain-24-cells` | `inflow_b24.bin` | `judge_approach_flow.py` |

Each run averages 1200 s after 600 s of spin-up. On one RTX 4090 the block runs took about 11, 51 and 146 minutes. The paths in the run settings are relative to `validation/`.

### Step 3: score and judge (CPU)

```
for n in 16 24 32; do
  uv run --with numpy python -I scripts/score_block_run.py runs/block-$n-cells inputs/block-$n-cells results/tunnel-data.json results/block-$n-cells-score.json
done
for n in 16 24; do
  uv run --with numpy python -I scripts/judge_approach_flow.py runs/empty-domain-$n-cells inputs/empty-domain-$n-cells results/tunnel-data.json results/empty-domain-$n-cells-approach.json
done
```

Both scripts take the run directory, the geometry folder, `results/tunnel-data.json` and an output file.

### How the calibration file was made

`inflow-calibration.json` was found in empty-domain runs at 16 cells and is used as given on all three grids. Its multipliers of the three standard deviations came from empty-domain runs judged with `judge_approach_flow.py`, which writes the next multipliers from the ones used and the ratio of the sheet value to the simulated value at the block position. Its multipliers of the mean speed came from two more such steps with the option `--r 1`. The calibration runs are not part of this recipe.

The second check (AIJ case D) is reported below but cannot be repeated from this repository: its geometry builder, run settings and scorer are not included.

## Results: case A

The earlier version scored 0.57 against 0.66. This version added a turbulent inflow calibrated on the tunnel's inflow sheet in empty-domain runs, the wall-stress ground and a Smagorinsky constant of 0.12.

| Case A, horizontal speed unless noted | Earlier version, 16 cells | 16 cells | 24 cells | 32 cells |
|---|---|---|---|---|
| B1: hit rate, 181 locations [95% interval] | 0.57 | 0.81 [0.74, 0.86] | 0.92 [0.87, 0.95] | 0.94 [0.89, 0.97] |
| B2: hit rate at pedestrian level | 0.50 | 0.83 | 0.93 | 0.95 |
| B3: hit rate behind the block, signed streamwise component | 0.35 | 0.68 | 0.85 | 0.94 |
| B4: hit rate of the weakest region | 0.39 (wake) | 0.70 (wake below roof) | 0.875 (sides) | 0.925 (sides) |
| B5: reattachment length at pedestrian level (tunnel 1.32 b) | 4.14 b | 1.67 b | 1.44 b | 1.32 b |
| Vertical component: hit rate, 181 locations | | 0.84 | 0.87 | 0.87 |

- Rows B1 to B4 are met on all three grids of this version. Row B5 is met at 24 and 32 cells and missed at 16 cells by 0.016 b. Row C1 is met on the pair 24 and 32 cells (the reattachment length changes by 8.3%) and not on the pair 16 and 24 (13.7%).
- The tunnel's reattachment length of 1.32 b is our linear interpolation between two measuring points 0.75 b apart; the papers cited for this case give 1.42 b for the same data.
- Row B6 was written for the steady flow of the earlier version and is not met as written on any grid (mean difference 0.040, 0.033 and 0.042 of the approach speed against 0.01, and largest difference 0.136, 0.136 and 0.151 against 0.05). After the results were known it was replaced by the standard error of the mean, at most 0.025 of the approach speed on average; by that rule it is met on all three grids (0.018, 0.018, 0.019).
- The three grids are not a pure refinement of one discrete problem: the shortest wavelength of the inflow, the floor of its length scale and the drag coefficient of the first cell above the ground change with the cell size.
- The approach flow was calibrated at 16 cells in an empty domain; that calibration passed on whole-window values with margins of about one standard error on the profile rows. At 32 cells no empty-domain run was made: the approach flow estimated from the block run was 5.1 to 5.3% high at 4 m above the ground against a mark of 5%.
- Read as a whole: the single block passes from 24 cells across it; at 16 cells four of five rows are met.

## Second case: AIJ case D, not tuned

AIJ case D (a 100 m tower among 78 city blocks, wind 0 degrees, 33 independent locations at 2 m above the ground) was declared before its run and scored once on a 2 m grid, with the binary and settings of case A and one calibration run of its own approach flow: hit rate 0.58 [0.41, 0.73] against 0.66, not met. The simulated wind in the streets is too slow (mean bias -0.085 of the approach speed, -0.144 in the 10 m streets).

A diagnostic afterwards, in a reduced domain at 1.25 m cells and not a scoring, gave 0.79 [0.62, 0.89]; the same reduced domain at 2 m gave 0.64. That diagnostic had no empty-domain run of its own, was normalised by the tabulated approach speed, ran in a reduced domain over a shorter window, and still left a bias of -0.091 in the 10 m streets.

The geometry of case D was read from the workbook's figure, because the AIJ geometry files could not be obtained.

## Limits

- Calibration and scoring on the same case. The inflow was calibrated in empty-domain runs of case A, judged at the block position and never on the points around the block; the correction of the mean inlet speed that closed that calibration was written after the declared adjustments had failed. The Smagorinsky constant 0.12 was chosen after a diagnostic run with the block, from reference studies; the controls show it was not the most favourable value (0.10 scored better on every row at 16 cells). Case A is therefore not a blind test. Case D was declared before its run and scored once, but its regions were chosen with the measured values in view; it did not pass.
- One pass rule (B6) was changed after the results were known.
- Reynolds-number independence is assumed, not shown: the lattice Reynolds number is about 1700 times the tunnel's.
- Mean wind only. The pass marks score the mean wind. Turbulence was compared without a pass mark (behind the block sigma_u is within 1% of the tunnel's at 16 cells; sigma_v and sigma_w are 8 to 17% above the tunnel's at 16 and 24 cells). Nothing here supports statements about gusts.
- The inflow carries no Reynolds shear stress, and its profile still changes along the approach (4.4% over the last two building heights at 16 cells).
- The side and top faces hold a fixed profile instead of the slip faces the AIJ guideline asks for. Building walls and roofs have no wall model.
- One wind direction; the hit-rate intervals treat neighbouring points as independent.
- Controls at 16 cells: 32-bit storage of the distributions changed B1 from 0.81 to 0.83 and the reattachment length from 1.67 to 1.55 b. The two-relaxation-time operator was not tried.

What this supports: the flow at the scale of a building on grids of 24 cells across it or finer. What it does not support: wind near the ground in streets and gaps of fewer than about ten to fifteen cells across.

The limits of the setup itself (boundaries, solids, inflow table, rendered images) are in [USAGE.md](USAGE.md#limits-of-the-setup).

## References

- Tominaga, Y., Mochida, A., Yoshie, R., Kataoka, H., Nozu, T., Yoshikawa, M., Shirasawa, T. (2008) AIJ guidelines for practical applications of CFD to pedestrian wind environment around buildings, J. Wind Eng. Ind. Aerodyn. 96, 1749-1761.
- Meng, Y., Hibi, K. (1998) Turbulent measurements of the flow field around a high-rise building, J. Wind Eng. Japan 76, 55-64.
