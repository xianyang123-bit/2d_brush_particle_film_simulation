"""Direct in-plane core S(q) on periodic reciprocal vectors; production only.
Requires numpy, matplotlib, gsd. Usage: python core_structure_factor.py RUN_DIR
"""
import argparse
import csv
import json
import os
from contextlib import nullcontext
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS'] = '2'
os.environ['OMP_NUM_THREADS'] = '2'
os.environ['MKL_NUM_THREADS'] = '2'
import numpy as np
import gsd.hoomd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('run', type=Path)
p.add_argument('--qmax', type=float, default=6.0)
a = p.parse_args()
out = a.run / 'core_structure_factor'
out.mkdir(exist_ok=True)
stages = json.loads((a.run / 'stages.json').read_text())
stage = next(s for s in stages if s['stage'] == 'production')
positions, steps, boxes = [], [], []
with gsd.hoomd.open(str(a.run / 'trajectory.gsd'), 'r') as traj:
    for f in traj:
        step = int(f.configuration.step)
        if stage['start'] < step <= stage['end']:
            mask = f.particles.typeid == f.particles.types.index('C')
            positions.append(f.particles.position[mask, :2].astype(float))
            steps.append(step)
            boxes.append(np.array(f.configuration.box, dtype=float))
xy = np.array(positions)
boxes = np.array(boxes)
assert len(xy) > 1 and np.isfinite(xy).all()
assert np.allclose(boxes, boxes[0], rtol=0, atol=1e-6)
assert np.allclose(boxes[:, 3:], 0)
Lx, Ly = boxes[0, :2]
assert np.isclose(Lx, Ly)
nframes, N, _ = xy.shape
q0 = 2 * np.pi / Lx
m = int(np.ceil(a.qmax / q0))
nx = np.arange(-m, m + 1)
ny = nx.copy()
gx, gy = np.meshgrid(nx, ny, indexing='ij')
qm = q0 * np.sqrt(gx**2 + gy**2)
# S(-q)=S(q): retain one member of each inversion pair.
keep = ((gx > 0) | ((gx == 0) & (gy > 0))) & (qm <= a.qmax)
q = qm[keep]
qx, qy = (q0 * gx[keep]), (q0 * gy[keep])
def amplitudes(pos):
    ex = np.exp(1j * q0 * nx[:, None] * pos[None, :, 0])
    ey = np.exp(1j * q0 * pos[:, None, 1] * ny[None, :])
    return ex @ ey

values = []
with nullcontext():
    for i, pos in enumerate(xy):
        rho = amplitudes(pos)
        assert np.isclose(abs(rho[m, m])**2 / N, N)
        values.append((abs(rho[keep])**2) / N)
        if (i + 1) % 20 == 0:
            print(f'Computed {i+1}/{nframes} frames', flush=True)
    # Independent direct-sum check and invariance to a wrapped translation.
    picks = np.linspace(0, len(q)-1, 23, dtype=int)
    direct = abs(np.exp(1j*(xy[0,:,0,None]*qx[picks] +
                            xy[0,:,1,None]*qy[picks])).sum(axis=0))**2 / N
    np.testing.assert_allclose(np.array(values)[0, picks], direct, atol=1e-10)
    shifted = (xy[0] + [7.321, -11.123] + [Lx/2, Ly/2]) % [Lx, Ly] - [Lx/2, Ly/2]
    np.testing.assert_allclose(abs(amplitudes(shifted)[keep])**2/N, values[0], atol=1e-9)
values = np.array(values)
assert np.isfinite(values).all() and (values >= 0).all()
np.savez_compressed(out/'core_positions_production.npz', xy=xy, steps=steps, box=boxes[0])

# Bin width one fundamental reciprocal spacing; report actual mean q of vectors.
binid = np.floor(q/q0).astype(int)
ids = np.unique(binid)
radial = np.array([values[:, binid == b].mean(axis=1) for b in ids]).T
centers = np.array([q[binid == b].mean() for b in ids])
counts = np.array([(binid == b).sum() for b in ids])
mean = radial.mean(axis=0)
blocks = np.array([r.mean(axis=0) for r in np.array_split(radial, 5)])
spread = blocks.std(axis=0, ddof=1)
with (out/'core_structure_factor.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['q_sigma_inverse','S_q','block_mean_std','independent_q_vectors','frames'])
    w.writerows(zip(centers, mean, spread, counts, [nframes]*len(mean)))
np.savez_compressed(out/'structure_factor_frames.npz', q=centers, S_by_frame=radial,
                    steps=steps, block_means=blocks, independent_q_vectors=counts)

# Exact reciprocal shells retain low-q resolution without combining shell radii.
n2 = gx[keep]**2 + gy[keep]**2
shells = np.unique(n2[q <= 0.6])
low_q = q0*np.sqrt(shells)
low_values = np.array([values[:, n2 == s].mean(axis=1) for s in shells]).T
low_mean = low_values.mean(axis=0)
low_blocks = np.array([r.mean(axis=0) for r in np.array_split(low_values,5)])
low_spread = low_blocks.std(axis=0, ddof=1)
with (out/'core_structure_factor_lowq_shells.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['q_sigma_inverse','S_q','block_mean_std','independent_q_vectors'])
    w.writerows(zip(low_q, low_mean, low_spread, [(n2==s).sum() for s in shells]))

plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False})
fig, ax = plt.subplots(1,2,figsize=(12,4.6),gridspec_kw={'width_ratios':[1.6,1]})
for axis, xx, yy, ss in zip(ax, [centers,low_q], [mean,low_mean], [spread,low_spread]):
    axis.axhline(1,color='gray',ls='--',lw=1,label='Uncorrelated reference: S = 1')
    axis.fill_between(xx,np.maximum(0,yy-ss),yy+ss,color='#007c91',alpha=.2,
                      label='±1 SD of 5 time-block means')
    axis.plot(xx,yy,'o-',color='#007c91',lw=1.4,ms=2.5,label='Production average')
    axis.set(xlabel=r'$q\sigma$ (angular wave number)',ylabel=r'$S_{\mathrm{core}}(q)$',ylim=(0,None))
    axis.grid(alpha=.18)
ax[0].set(xlim=(0,a.qmax),title='In-plane core structure factor')
ax[1].set(xlim=(0,.6),ylim=(0,1.12*np.max(low_mean+low_spread)),title='Small q: exact reciprocal shells')
ax[0].legend(fontsize=8,loc='upper right')
fig.suptitle(f'{N} nanoparticle cores · {nframes} production frames · steps {steps[0]:,}–{steps[-1]:,}',fontsize=12)
fig.text(.5,.01,'q = 0 excluded. Shading shows time-block variability, not a confidence interval.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.94])
for ext in ['png','pdf']:
    fig.savefig(out/f'core_structure_factor.{ext}',dpi=220)
plt.close(fig)
peak = np.argmax(mean)
summary = dict(cores=N,frames=nframes,first_step=steps[0],last_step=steps[-1],Lx=Lx,
    q_min=q0,q_max=a.qmax,unique_reciprocal_vectors=len(q),radial_bin_width=q0,
    S_at_smallest_q=float(low_mean[0]),smallest_q_block_std=float(low_spread[0]),
    radial_peak_q=float(centers[peak]),radial_peak_S=float(mean[peak]),
    validation='Direct-sum agreement, periodic translation invariance, S(0)=N, finite and nonnegative checks passed.')
(out/'summary.json').write_text(json.dumps(summary,indent=2))
(out/'METHODS.md').write_text('''# Core structure factor

Computed from core type C only, projecting positions onto the periodic xy plane.
S(q) = < |sum_j exp(i q dot r_j)|^2 / N > over production frames.
q = 2 pi (nx/Lx, ny/Ly); q is the angular wave number (not cycles per length).
Only production frames strictly after step 2,600,000 are used. Each frame has equal weight.
q=0 is excluded. Inversion pairs q and -q have identical intensity and count once.
The full curve averages reciprocal vectors in bins of width 2 pi/Lx; x values are
the mean magnitude of included vectors. The small-q curve uses exact reciprocal shells.
No particle form factor or polymer contribution is included. This is a 2D in-plane
number-density structure factor of the finite-thickness film, not spherical 3D scattering.

Shading is one sample standard deviation among five consecutive 20-frame block means.
Blocks may remain correlated; this is descriptive variability, not a confidence interval.
The 400-core finite box limits small-q resolution. Completion, a low S(q), or this
curve alone cannot establish equilibration or hyperuniformity.

Definition reference: https://freud.readthedocs.io/en/latest/gettingstarted/examples/module_intros/diffraction.StaticStructureFactor.html
The calculation here uses a direct NumPy sum on the 2D reciprocal lattice.
''')
print(json.dumps(summary,indent=2))
