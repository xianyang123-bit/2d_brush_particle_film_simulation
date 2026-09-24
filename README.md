# Polymer-grafted nanoparticle film in HOOMD-blue

This repository contains a HOOMD-blue 7.1 model of polymer-grafted nanoparticle
cores confined between two repulsive walls. It was developed with the simulation
model in Chremos and Douglas, *Annalen der Physik* **529**, 1600342 (2017),
[doi:10.1002/andp.201600342](https://doi.org/10.1002/andp.201600342), as a
reference. The paper studies a bulk three-dimensional material; this code is a
finite-thickness film adaptation and is not an exact reproduction of its phase
diagram.

Each diameter-3 core has eight grafted chains. A chain contains one rigid surface
anchor and nine flexible beads. Core-core and core-polymer interactions are
repulsive expanded Lennard-Jones interactions. Polymer-polymer interactions use
an attractive Lennard-Jones potential with cutoff 2.5. The film is periodic in
the lateral directions and confined by walls separated by 4 reduced length units.

## Environment

HOOMD-blue should be run under Linux, WSL2, or macOS. Create the environment with:

```bash
conda env create -f environment.yml
conda activate gnp2d
```

The environment pins HOOMD-blue 7.1.x. Despite the historical environment name,
the current model is a three-dimensional slit film, not a strict 2D simulation.

## Run the 400-core film

`run_film400_long.sh` launches the tested CPU configuration in the background and
stores each run in a new directory. It expects the environment installed by the
original helper at `~/.local/share/gnp-hoomd/env`. To use an activated conda
environment instead, run the Python command directly:

```bash
python gnp_film.py \
  --nside 20 --gap 4 --kT 1.0 --dt 0.0005 \
  --tau-thermostat 0.2 --tau-barostat 50 \
  --relax-steps 100000 \
  --prep-pressure 0.1 --prep-steps 500000 \
  --pressure 0 --npt-steps 1000000 \
  --equil-steps 1000000 --production-steps 1000000 \
  --period 10000 --checkpoint-period 100000 \
  --out film400

python analyze_film.py film400
```

Every simulation requires a new output directory. This prevents accidental
overwriting of previous trajectories.

## Structure-factor analysis

The analysis uses nanoparticle core centers projected into the periodic xy plane:

```text
S(q) = < |sum_j exp(i q . r_j)|^2 / N_core >
q = 2 pi (nx/Lx, ny/Ly)
```

Run:

```bash
python core_structure_factor.py film400
python compare_structure_factor_frames.py film400/core_structure_factor
python plot_smallest_q_frames.py film400/core_structure_factor
```

The calculation selects production frames from `stages.json`, excludes q=0, and
averages reciprocal vectors radially. It also exports exact low-q shells and
comparisons of single-frame, 5-frame, and 20-frame averages.

In the completed 400-core run, the principal peak was at q sigma approximately
1.02 with S(q) approximately 2.06. Suppressed
small-q intensity alone does not establish hyperuniformity: finite-size scaling,
stationarity, equilibration, and disorder must also be checked.

## Files

- `gnp_film.py`: constructs and runs the confined film.
- `analyze_film.py`: final geometry and thermodynamic diagnostics.
- `core_structure_factor.py`: direct in-plane core structure factor.
- `compare_structure_factor_frames.py`: 1-, 5-, and 20-frame averaging comparison.
- `plot_smallest_q_frames.py`: smallest reciprocal-shell time series.
- `FILM_README.md`: detailed model assumptions and force definitions.
- `FILM400_RUN.md`: notes for the tested 400-core run.

Large GSD trajectories and checkpoints are intentionally excluded from version
control.
