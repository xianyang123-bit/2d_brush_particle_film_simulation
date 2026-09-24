"""Inspect completed wall-confined film runs, without claiming equilibration."""
import argparse
import csv
import json
from pathlib import Path
import gsd.hoomd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('run', type=Path)
a = p.parse_args()
config = json.loads((a.run/'parameters.json').read_text())
with gsd.hoomd.open(str(a.run/'production_end.gsd'), 'r') as traj:
    frame = traj[-1]
pos = frame.particles.position
tid = frame.particles.typeid
core = tid == frame.particles.types.index('C')
H = config['gap']
box = np.asarray(frame.configuration.box[:3])
dv = pos[frame.bonds.group[:,0]] - pos[frame.bonds.group[:,1]]
dv -= box*np.rint(dv/box)
bond_lengths = np.linalg.norm(dv,axis=1)
with (a.run/'thermo.csv').open() as f:
    rows = list(csv.DictReader(f))
numeric = lambda key: np.array([float(r[key]) for r in rows])
summary = dict(hoomd_version=(a.run/'hoomd_version.txt').read_text().strip(),
    dimensions=int(frame.configuration.dimensions), cores=int(core.sum()), sites=len(pos),
    wall_separation=H, steps=int(frame.configuration.step),
    physical_time=float(frame.configuration.step)*config['dt'],
    final_core_z_min=float(pos[core,2].min()), final_core_z_max=float(pos[core,2].max()),
    final_bead_z_min=float(pos[~core,2].min()), final_bead_z_max=float(pos[~core,2].max()),
    final_bond_min=float(bond_lengths.min()), final_bond_max=float(bond_lengths.max()),
    all_logged_values_finite=bool(all(np.isfinite(numeric(key)).all() for key in rows[0] if key!='stage')),
    no_logged_wall_crossing=bool(max(numeric('max_abs_z_core').max(),numeric('max_abs_z_bead').max()) < H/2),
    final_temperature=float(numeric('kT')[-1]), final_projected_core_area_fraction=float(numeric('phi_core_projected')[-1]),
    interpretation='Diagnostic run only; equilibration and a homogeneous monolayer are not established.')
(a.run/'run_check.json').write_text(json.dumps(summary,indent=2))
colors = np.array(['#244a74','#ed9c40','#4eada7'])[tid]
sizes = np.where(core,config['diameter'],1.)
fig,(top,side)=plt.subplots(2,1,figsize=(10,10),gridspec_kw={'height_ratios':[3,1]})
for ax,axes in [(top,(0,1)),(side,(0,2))]:
    ax.add_collection(PatchCollection([Circle((q[axes[0]],q[axes[1]]),s/2) for q,s in zip(pos,sizes)],
                                     facecolors=colors,edgecolors='none',alpha=.8))
    ax.set_xlim(-box[0]/2,box[0]/2)
top.set(ylim=(-box[1]/2,box[1]/2),aspect='equal',xlabel='x / sigma',ylabel='y / sigma',
        title=f'Film diagnostic run: {core.sum()} cores, step {frame.configuration.step}\nTop view; not asserted to be equilibrated')
side.set(ylim=(-H/2-1,H/2+1),xlabel='x / sigma',ylabel='z / sigma',title='Side projection (vertical scale expanded)')
for z in [-H/2,H/2]: side.axhline(z,color='#333333',lw=2)
fig.tight_layout();fig.savefig(a.run/'final_film.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(10,7))
for ax,key,label in zip(axs.flat,['kT','P_parallel_slit','phi_core_projected','max_bond_error'],
                       ['Temperature','Lateral pressure (slit volume)','Projected core area fraction','Largest bond length deviation']):
    ax.plot(numeric('time'),numeric(key),lw=1)
    ax.set(xlabel='time / tau',ylabel=label)
    ax.grid(alpha=.2)
fig.suptitle('Diagnostics include preparation and production')
fig.tight_layout();fig.savefig(a.run/'diagnostics.png',dpi=150);plt.close(fig)
print(json.dumps(summary,indent=2))
