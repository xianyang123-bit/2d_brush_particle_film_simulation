# Finite-thickness film between two walls

Use `gnp_film.py` for the clarified request. The earlier `gnp2d.py` describes a different geometry: a mathematical plane in which every center has z = 0.

The film has a **3D simulation box**, periodic boundaries in x and y, and two repulsive planes at z = -H/2 and +H/2. All cores, grafts and free chain beads feel the walls. Positions have three translational coordinates; cores rotate about all three axes. The default H = 4 sigma with core diameter D = 3 sigma is a chosen monolayer slit, not a value from the paper. A monolayer describes the nanoparticle cores; the polymer beads occupy the finite slit thickness.

Run in an existing x86_64 Ubuntu/WSL terminal:

```bash
bash '/mnt/c/Users/xianyang chen/Documents/Codex/2026-09-22/h/outputs/install_and_run_film.sh'
```

The script installs a CPU HOOMD 7.1.x environment under `${XDG_DATA_HOME:-$HOME/.local/share}/gnp-hoomd/env`, using conda-forge via micromamba. It does not require sudo or modify your shell profile. It first executes a 4-core, 1000-step smoke test; then a 16-core, 40000-step diagnostic run. Each launch creates a fresh `film_run_XXXXXXXX` output directory, with installation log and status. `latest_film_run.txt` identifies it. The script stops if a command fails. Initial package downloads require internet access.

The current Codex process cannot launch WSL: Windows returns `Wsl/EnumerateDistros/Service/E_ACCESSDENIED`. Consequently this installer must be launched from the user's working Ubuntu terminal. Prepared input files and previews alone are not evidence that HOOMD was installed or that dynamics ran. Consult the actual run status and `run_check.json` when available.

## Wall and graft choices

Each wall uses repulsive LJ truncated at its minimum, with epsilon_wall = 1. The wall sigma is D/2 for core centers, and 1/2 for polymer bead centers. Thus wall forces start at center-wall distances 1.68369 and 0.561231 sigma, respectively. With H = 4, the zero-wall-force band for core centers is about -0.3163 < z < 0.3163. These are soft potentials, not mathematically hard boundaries: the code checks sampled states for escaped particles, and a timestep convergence test remains necessary.

The 8 grafts lie on the core surface at z = +/- D/(2*sqrt(3)), alternating around the sphere; their xy radius is (D/2)*sqrt(2/3). This symmetric 3D arrangement has zero center-of-mass offset and diagonal isotropic inertia. It is a chosen graft pattern, not an extraction of the paper's exact graft coordinates. Initially the flexible parts extend horizontally at the graft heights so they fit inside the slit. Three-dimensional velocities and rotations allow them to explore the slit during dynamics.

The bare core mass remains 27, with rigid-body mass 35 including its graft beads. With D = 3, the solid-sphere core plus point grafts gives Ixx = Iyy = Izz = 36.3. Particle interactions and chain parameters otherwise match the earlier documented adaptation: 8 chains of 10 beads including the graft, energy-shifted attractive bead LJ at cutoff 2.5, harmonic k_HOOMD = 10000, and expanded WCA core forces. **The paper's ambiguous core-bead shift remains unresolved**: this implementation chooses delta_nb = (D - 1)/2 = 1, while page 3 prints D - 1/2 = 2.5. See `README.md` for the distinction and bonded-pair exclusion choice.

The numerical box is padded outside both walls so periodic replicas cannot interact across z: Lz - H exceeds the maximum pair cutoff. Provided particles remain inside the walls, this produces an isolated slit even though HOOMD uses a periodic box internally. The outside region is empty padding, not accessible film volume. The wall separation and Lz stay fixed when the lateral box area changes.

## Pressure and run interpretation

Only the x and y lengths are barostatted, coupled equally. HOOMD's stress tensor uses the full padded volume Lx*Ly*Lz, so a desired slit lateral pressure P is supplied to the barostat as P*H/Lz. The CSV reports `(Pxx + Pyy)/2 * Lz/H`, normalized by the physical slit volume Lx*Ly*H. It is not the isotropic pressure or the wall-normal pressure. Units are epsilon/sigma^3. The target is zero by default after a positive-pressure preparation stage. Fixed-area production follows.

The automated diagnostic is deliberately short: 5000 fixed-area relaxation steps, 20000 positive-pressure preparation steps, 5000 target-pressure steps, 5000 fixed-area equilibration steps, and 5000 production steps, all at dt = 0.0005. **It does not establish equilibrium, a homogeneous monolayer, or hyperuniformity.** The paper studied bulk 3D materials; its phase diagram and zero-pressure state need not survive confinement. Inspect the time traces, z distributions, density, and chain relaxation, and increase durations as required.

For larger runs after the small test succeeds, use the installed Linux interpreter, e.g.:

```bash
"${XDG_DATA_HOME:-$HOME/.local/share}/gnp-hoomd/env/bin/python" gnp_film.py --nside 20 --gap 4 --out film400
```

This produces 400 cores and 32000 polymer beads. `--gap` changes the wall separation; larger values can permit multiple core layers. Use the shortest timestep that passes convergence testing rather than assuming the paper's 0.005 is suitable with added walls.

After the automated run, `film16/final_film.png` shows top and side views, `diagnostics.png` shows time traces, and `run_check.json` records dimensions, particle counts, wall containment at logged times, bond lengths and thermodynamics. Full trajectory and stage endpoints are saved as GSD. The side view expands the vertical scale to show the thin slit clearly.

References: [HOOMD wall LJ](https://hoomd-blue.readthedocs.io/en/stable/hoomd/md/external/wall/lj.html), [plane geometry](https://hoomd-blue.readthedocs.io/en/stable/hoomd/wall/plane.html), [micromamba installation](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).
