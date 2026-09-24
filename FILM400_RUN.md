# Longer 400-core film run

Launch `run_film400_long.sh` from Ubuntu. It reuses the installed CPU HOOMD environment, copies the simulation and analysis source into a new run directory, and launches a background process. Closing the terminal should not terminate it, but shutting down WSL, rebooting, or suspending the computer prevents continuous progress. Running the launcher again starts a separate new job.

The model remains a 3D slit film with H = 4 sigma and D = 3 sigma, 8 grafted chains per core and 10 beads per chain. There are 400 nanoparticle cores, 32000 polymer beads, and 32400 interaction sites. Force-field assumptions are unchanged from FILM_README.md.

| Stage | Steps | Reduced time |
|---|---:|---:|
| Fixed-area relaxation | 100000 | 50 |
| Preparation at lateral P = 0.1 | 500000 | 250 |
| Lateral P = 0 equilibration | 1000000 | 500 |
| Fixed-area equilibration | 1000000 | 500 |
| Fixed-area production | 1000000 | 500 |
| Total | 3600000 | 1800 |

The timestep remains 0.0005. This is 25 times as many nanoparticle cores and 90 times the physical duration of the completed 16-core diagnostic. CPU wall time is not yet measured for this configuration; expect a substantially longer job.

To moderate the sharp compression and temperature oscillations seen in the short run, preparation pressure is lowered from 0.5 to 0.1, barostat coupling time is increased from 10 to 50, and thermostat coupling time is decreased from 1 to 0.2. These are trial numerical choices, not established optimal settings. Temperature, density and structural relaxation still need inspection; completing the longer schedule does not automatically establish equilibration.

Trajectories and CSV diagnostics are sampled every 10000 steps. Every 100000 steps the code flushes the trajectory, saves `checkpoint.gsd`, writes `progress.json` and prints progress to `run.log`. The checkpoint stores particle state, not full thermostat/barostat state, so an exact continuation is not implemented. `latest_film400_run.txt` identifies the output directory; `pid.txt` holds its Linux process ID. At completion the analysis produces final views, diagnostics and `run_check.json`.
