"""3D polymer-grafted nanoparticle monolayer between two repulsive walls.

Targets HOOMD-blue 7.x. --init-only needs only NumPy and Matplotlib.
"""
import argparse
import csv
import itertools
import json
import time
from pathlib import Path

import numpy as np


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nside', type=int, default=4, help='N cores = nside squared; use 20 for 400')
    p.add_argument('--chains', type=int, default=8)
    p.add_argument('--beads', type=int, default=10, help='Includes the rigid first bead')
    p.add_argument('--diameter', type=float, default=3.0)
    p.add_argument('--gap', type=float, default=4.0, help='Physical wall separation in bead sigma units')
    p.add_argument('--graft-radius', type=float, default=None)
    p.add_argument('--delta-nb', type=float, default=None)
    p.add_argument('--kT', type=float, default=1.0)
    p.add_argument('--dt', type=float, default=0.0005)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--relax-steps', type=int, default=20000)
    p.add_argument('--prep-steps', type=int, default=200000)
    p.add_argument('--prep-pressure', type=float, default=0.5)
    p.add_argument('--npt-steps', type=int, default=200000)
    p.add_argument('--pressure', type=float, default=0.0)
    p.add_argument('--equil-steps', type=int, default=100000)
    p.add_argument('--production-steps', type=int, default=100000)
    p.add_argument('--period', type=int, default=1000)
    p.add_argument('--checkpoint-period', type=int, default=100000)
    p.add_argument('--tau-thermostat', type=float, default=1.0)
    p.add_argument('--tau-barostat', type=float, default=10.0)
    p.add_argument('--exclude-bonded-lj', action='store_true')
    p.add_argument('--device', choices=['cpu', 'gpu'], default='cpu')
    p.add_argument('--init-only', action='store_true')
    p.add_argument('--out', type=Path, default=Path('run_film'))
    a = p.parse_args()
    if a.nside < 1 or a.chains < 2 or a.beads < 2 or a.diameter < 1:
        p.error('Require nside >= 1, chains >= 2, beads >= 2, diameter >= 1.')
    if min(a.kT, a.dt, a.tau_thermostat, a.tau_barostat) <= 0 or min(a.period, a.checkpoint_period) < 1:
        p.error('Temperature, timestep, coupling times, and period must be positive.')
    if any(getattr(a, s) < 0 for s in ['relax_steps', 'prep_steps', 'npt_steps', 'equil_steps', 'production_steps']):
        p.error('Step counts cannot be negative.')
    if not 0 <= a.seed <= 65535:
        p.error('Use a seed between 0 and 65535.')
    if a.graft_radius is None:
        a.graft_radius = a.diameter / 2
    if a.delta_nb is None:
        a.delta_nb = (a.diameter - 1) / 2
    if a.graft_radius <= 0 or a.delta_nb < 0:
        p.error('Require graft-radius > 0 and delta-nb >= 0.')
    if a.chains != 8:
        p.error('The film initializer currently supports exactly 8 symmetric surface grafts.')
    if a.gap <= max(2**(1/6)*a.diameter, 2*(a.graft_radius/np.sqrt(3)+2**(1/6)/2)):
        p.error('Gap is too narrow for the initial cores and graft beads without wall overlap.')
    first_flexible_radius = np.sqrt((a.graft_radius*np.sqrt(2/3)+1)**2 + a.graft_radius**2/3)
    if first_flexible_radius - a.delta_nb < 1:
        p.error('First flexible bead is inside the core repulsion: increase graft-radius or reduce delta-nb.')
    return a


def build(a):
    """Symmetric 3D grafts; chains initially extend horizontally inside the slit."""
    nc = a.nside ** 2
    n = nc * (1 + a.chains * a.beads)
    reach = a.graft_radius + a.beads - 1
    spacing = max(2 * reach + 4, 2 * (a.diameter - 1 + 2**(1/6) + 0.4) + 1)
    L = a.nside * spacing
    theta = np.arange(a.chains) * 2 * np.pi / a.chains
    rho = a.graft_radius * np.sqrt(2/3)
    heights = np.where(np.arange(a.chains) % 2 == 0, 1., -1.) * a.graft_radius/np.sqrt(3)
    local = np.column_stack((rho*np.cos(theta), rho*np.sin(theta), heights))
    pos = np.zeros((n, 3))
    orientation = np.tile([1., 0., 0., 0.], (n, 1))
    tid = np.full(n, 2, dtype=np.int32)  # C = core, A = anchor, B = flexible bead
    body = np.full(n, -1, dtype=np.int32)
    mass = np.ones(n)
    inertia = np.zeros((n, 3))
    diameter = np.ones(n)
    bonds = []
    rng = np.random.default_rng(a.seed)
    tid[:nc] = 0
    body[:nc] = np.arange(nc)
    # Central mass represents core + rigid anchors; anchor masses are not integrated.
    mass[:nc] = a.diameter**3 + a.chains
    # Solid spherical core + symmetric point anchors: three active rotation axes.
    inertia[:nc, :] = (2/5)*a.diameter**3*(a.diameter/2)**2 + (2/3)*a.chains*a.graft_radius**2
    diameter[:nc] = a.diameter
    cursor = nc
    for c, (ix, iy) in enumerate(itertools.product(range(a.nside), repeat=2)):
        pos[c] = [(ix + .5) * spacing - L/2, (iy + .5) * spacing - L/2, 0]
        angle = rng.uniform(0, 2*np.pi)
        orientation[c] = [np.cos(angle/2), 0, 0, np.sin(angle/2)]
        for j in range(a.chains):
            direction = np.array([np.cos(theta[j]+angle), np.sin(theta[j]+angle), 0])
            anchor = cursor
            for b in range(a.beads):
                pos[cursor] = pos[c] + (rho + b) * direction + [0, 0, heights[j]]
                if b == 0:
                    tid[cursor] = 1
                    body[cursor] = c
                    orientation[cursor] = orientation[c]
                else:
                    bonds.append((cursor-1, cursor))
                cursor += 1
            assert body[anchor] == c
    result = dict(position=pos, orientation=orientation, typeid=tid, body=body,
                  mass=mass, moment_inertia=inertia, diameter=diameter,
                  bonds=np.asarray(bonds, dtype=np.int32), local=local, L=np.array(L),
                  Lz=np.array(a.gap + 2*max(2.5, a.diameter-1+2**(1/6), a.delta_nb+2**(1/6)) + 2))
    validate(a, result)
    return result


def validate(a, d):
    nc = a.nside**2
    assert len(d['position']) == nc * (1 + a.chains*a.beads)
    assert len(d['bonds']) == nc*a.chains*(a.beads-1)
    assert np.max(np.abs(d['position'][:, 2])) < a.gap/2
    assert np.all(d['moment_inertia'][:nc] > 0)
    assert np.max(np.abs(d['position'][:, :2])) < float(d['L'])/2
    lengths = np.linalg.norm(d['position'][d['bonds'][:, 0]] - d['position'][d['bonds'][:, 1]], axis=1)
    assert np.allclose(lengths, 1)
    for c in range(nc):
        anchors = np.flatnonzero((d['body'] == c) & (d['typeid'] == 1))
        assert len(anchors) == a.chains and np.all(anchors > c)
        assert np.allclose(np.linalg.norm(d['position'][anchors]-d['position'][c], axis=1), a.graft_radius)
    # All flexible beads have independent translational DOF.
    assert np.all(d['body'][d['typeid'] == 2] == -1)


def preview(a, d):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import LineCollection, PatchCollection
    fig, (ax, side) = plt.subplots(2, 1, figsize=(9, 9), gridspec_kw={'height_ratios': [3, 1]})
    ax.add_collection(LineCollection(d['position'][d['bonds'], :2], colors='#537f8f', linewidths=.6))
    colors = ['#244a74', '#ed9c40', '#4eada7']
    patches = [Circle(p[:2], diam/2) for p, diam in zip(d['position'], d['diameter'])]
    ax.add_collection(PatchCollection(patches, facecolors=[colors[t] for t in d['typeid']], edgecolors='none'))
    L = float(d['L'])
    ax.set(xlim=(-L/2, L/2), ylim=(-L/2, L/2), aspect='equal', xlabel='x / sigma', ylabel='y / sigma',
           title=f'Wall-confined film: {a.nside**2} nanoparticles (initial geometry)\n'
                 f'Top view; {a.chains} chains/core, {a.beads} beads/chain')
    # A single row of cores avoids overplotting in the side projection.
    select = d['position'][:, 1] < -L/2 + (L/a.nside)
    for t in range(3):
        ids = select & (d['typeid'] == t)
        side.add_collection(PatchCollection([Circle((p[0], p[2]), diam/2)
            for p, diam in zip(d['position'][ids], d['diameter'][ids])], facecolor=colors[t], edgecolors='none'))
    for z in [-a.gap/2, a.gap/2]:
        side.axhline(z, color='#333333', linewidth=2)
    side.set(xlim=(-L/2,L/2), ylim=(-a.gap/2-1, a.gap/2+1),
             xlabel='x / sigma', ylabel='z / sigma', title=f'Side view of first row (vertical scale expanded); H = {a.gap:g} sigma')
    fig.tight_layout()
    fig.savefig(a.out/'initial.png', dpi=170)
    plt.close(fig)


def run(a, d):
    import hoomd
    device = hoomd.device.GPU() if a.device == 'gpu' else hoomd.device.CPU()
    if device.communicator.num_ranks != 1:
        raise RuntimeError('This starter uses one MPI rank; run directly without mpirun.')
    sim = hoomd.Simulation(device=device, seed=a.seed)
    (a.out/'hoomd_version.txt').write_text(hoomd.version.version + '\n')
    snap = hoomd.Snapshot()
    L = float(d['L'])
    Lz = float(d['Lz'])
    snap.configuration.box = [L, L, Lz, 0, 0, 0]
    snap.particles.N = len(d['position'])
    snap.particles.types = ['C', 'A', 'B']
    for key in ['position', 'orientation', 'typeid', 'body', 'mass', 'moment_inertia', 'diameter']:
        getattr(snap.particles, key)[:] = d[key]
    snap.bonds.N = len(d['bonds'])
    snap.bonds.types = ['polymer']
    snap.bonds.typeid[:] = 0
    snap.bonds.group[:] = d['bonds']
    sim.create_state_from_snapshot(snap)
    rigid = hoomd.md.constrain.Rigid()
    rigid.body['C'] = dict(constituent_types=['A']*a.chains,
                           positions=d['local'].tolist(),
                           orientations=[(1., 0., 0., 0.)]*a.chains)
    moving = hoomd.filter.Rigid(('center', 'free'))
    exclusions = ['body'] + (['bond'] if a.exclude_bonded_lj else [])
    nl = hoomd.md.nlist.Cell(buffer=0.4, exclusions=exclusions)
    lj = hoomd.md.pair.LJ(nlist=nl, default_r_cut=0., mode='shift')
    core = hoomd.md.pair.ExpandedLJ(nlist=nl, default_r_cut=0., mode='shift')
    for pair in itertools.combinations_with_replacement(['C', 'A', 'B'], 2):
        lj.params[pair] = dict(epsilon=0., sigma=1.)
        core.params[pair] = dict(epsilon=0., sigma=1., delta=0.)
        if 'C' in pair:
            delta = a.diameter - 1 if pair == ('C', 'C') else a.delta_nb
            core.params[pair] = dict(epsilon=1., sigma=1., delta=delta)
            core.r_cut[pair] = delta + 2**(1/6)
        else:
            lj.params[pair] = dict(epsilon=1., sigma=1.)
            lj.r_cut[pair] = 2.5
    spring = hoomd.md.bond.Harmonic()
    spring.params['polymer'] = dict(k=10000., r0=1.)
    wall = hoomd.md.external.wall.LJ(walls=[
        hoomd.wall.Plane(origin=(0,0,-a.gap/2), normal=(0,0,1)),
        hoomd.wall.Plane(origin=(0,0,a.gap/2), normal=(0,0,-1))])
    for t in ['C', 'A', 'B']:
        radius = a.diameter/2 if t == 'C' else .5
        wall.params[t] = dict(epsilon=1., sigma=radius,
                              r_cut=2**(1/6)*radius, r_extrap=0.)

    def nvt():
        return hoomd.md.methods.ConstantVolume(filter=moving,
            thermostat=hoomd.md.methods.thermostats.MTTK(kT=a.kT, tau=a.tau_thermostat))

    integrator = hoomd.md.Integrator(dt=a.dt, integrate_rotational_dof=True,
        rigid=rigid, forces=[lj, core, spring, wall], methods=[nvt()])
    sim.operations.integrator = integrator
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=moving)
    sim.operations.computes.append(thermo)
    sim.always_compute_pressure = True
    sim.state.thermalize_particle_momenta(filter=moving, kT=a.kT)
    sim.run(0)
    if sim.state.box.dimensions != 3:
        raise RuntimeError('A finite-thickness film requires a 3D box.')
    log = hoomd.logging.Logger(categories=['scalar'])
    log.add(thermo, quantities=['kinetic_temperature', 'potential_energy', 'pressure'])
    trajectory = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(a.period),
        filename=str(a.out/'trajectory.gsd'), mode='xb', dynamic=['property', 'momentum'], logger=log)
    trajectory.write_diameter = True
    sim.operations.writers.append(trajectory)
    hoomd.write.GSD.write(state=sim.state, filename=str(a.out/'initial.gsd'), mode='xb')
    stages = []
    with (a.out/'thermo.csv').open('x', newline='') as fp:
        writer = csv.writer(fp)
        writer.writerow(['stage', 'step', 'time', 'kT', 'U', 'P_parallel_slit', 'Lx', 'Ly',
                         'phi_core_projected', 'max_abs_z_core', 'max_abs_z_bead', 'max_bond_error'])

        def stage(name, steps):
            start = sim.timestep
            stop = start + steps
            wall_start = time.monotonic()
            next_checkpoint = start + a.checkpoint_period
            while sim.timestep < stop:
                sim.run(min(a.period, stop - sim.timestep))
                box = sim.state.box
                current = sim.state.get_snapshot()
                pxyz = current.particles.position
                types = current.particles.typeid
                zcore = np.max(np.abs(pxyz[types == 0, 2]))
                zbead = np.max(np.abs(pxyz[types != 0, 2]))
                if max(zcore, zbead) >= a.gap/2:
                    raise RuntimeError('A particle crossed a wall. Stop and reduce timestep; this is invalid data.')
                v = pxyz[d['bonds'][:, 0]] - pxyz[d['bonds'][:, 1]]
                lengths = np.array([box.Lx, box.Ly, box.Lz])
                v -= lengths*np.rint(v/lengths)
                bond_error = np.max(np.abs(np.linalg.norm(v, axis=1)-1))
                # HOOMD pressure uses padded box volume; report physical slit volume.
                tensor = thermo.pressure_tensor
                lateral = .5*(tensor[0]+tensor[3])*Lz/a.gap
                row = [name, sim.timestep, sim.timestep*a.dt, thermo.kinetic_temperature,
                       thermo.potential_energy, lateral, box.Lx, box.Ly,
                       a.nside**2*np.pi*(a.diameter/2)**2/(box.Lx*box.Ly), zcore, zbead, bond_error]
                if not np.isfinite(row[1:]).all():
                    raise RuntimeError('Nonfinite thermodynamics: reduce timestep and inspect configuration.')
                writer.writerow(row)
                fp.flush()
                if sim.timestep >= next_checkpoint or sim.timestep == stop:
                    trajectory.flush()
                    # Snapshot checkpoint, not an exact thermostat/barostat restart.
                    temporary = a.out/'checkpoint.tmp.gsd'
                    hoomd.write.GSD.write(state=sim.state, filename=str(temporary), mode='wb')
                    temporary.replace(a.out/'checkpoint.gsd')
                    elapsed = time.monotonic()-wall_start
                    rate = (sim.timestep-start)/max(elapsed,1e-9)
                    progress = dict(stage=name, step=sim.timestep, stage_end=stop,
                        steps_per_second=rate, stage_eta_seconds=(stop-sim.timestep)/max(rate,1e-9),
                        temperature=float(thermo.kinetic_temperature), Lx=box.Lx)
                    (a.out/'progress.json').write_text(json.dumps(progress,indent=2))
                    print(f'{name}: {sim.timestep}/{stop}, T={thermo.kinetic_temperature:.4f}, '
                          f'L={box.Lx:.4f}, {rate:.1f} steps/s', flush=True)
                    next_checkpoint = sim.timestep + a.checkpoint_period
            trajectory.flush()
            hoomd.write.GSD.write(state=sim.state, filename=str(a.out/f'{name}_end.gsd'), mode='xb')
            stages.append(dict(stage=name, start=start, end=sim.timestep))
            (a.out/'stages.json').write_text(json.dumps(stages, indent=2))
            print(f'{name}: step={sim.timestep}, L={sim.state.box.Lx:.5g}', flush=True)

        stage('relax', a.relax_steps)
        if a.nside > 1:
            for name, steps, pressure in [('preparation', a.prep_steps, a.prep_pressure),
                                          ('npt', a.npt_steps, a.pressure)]:
                if steps:
                    integrator.methods = [hoomd.md.methods.ConstantPressure(
                        filter=moving, S=pressure*a.gap/Lz, tauS=a.tau_barostat, couple='xy',
                        box_dof=[True, True, False, False, False, False],
                        thermostat=hoomd.md.methods.thermostats.MTTK(kT=a.kT, tau=a.tau_thermostat))]
                    stage(name, steps)
        integrator.methods = [nvt()]
        stage('nvt_equil', a.equil_steps)
        stage('production', a.production_steps)
    (a.out/'hoomd_version.txt').write_text(hoomd.version.version + '\n')
    print('Finished requested steps. Assess equilibration before interpreting production data.')


def main():
    a = arguments()
    d = build(a)
    # Fail rather than overwrite a previous simulation.
    a.out.mkdir(parents=True, exist_ok=False)
    metadata = dict(vars(a), out=str(a.out), source_doi='10.1002/andp.201600342',
                    dimensions=3, geometry='slit film', harmonic_k_hoomd=10000, core_mass=a.diameter**3,
                    rigid_mass=a.diameter**3+a.chains, initial_L=float(d['L']), Lz=float(d['Lz']))
    (a.out/'parameters.json').write_text(json.dumps(metadata, indent=2))
    np.savez_compressed(a.out/'initial.npz', **d)
    preview(a, d)
    print(f"Built {a.nside**2} cores, {len(d['position'])} sites, {len(d['bonds'])} bonds.")
    if not a.init_only:
        run(a, d)


if __name__ == '__main__':
    main()
