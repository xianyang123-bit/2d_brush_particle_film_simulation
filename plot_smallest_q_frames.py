"""Plot exact smallest reciprocal-shell core structure factor by production frame."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('analysis_dir', type=Path)
a = p.parse_args()
d = np.load(a.analysis_dir / 'core_positions_production.npz')
xy, steps, box = d['xy'], d['steps'], d['box']
assert np.isclose(box[0], box[1])
qmin = 2*np.pi/box[0]
# Two independent directions. Negative vectors have identical intensities.
directional = np.abs(np.exp(1j*qmin*xy).sum(axis=1))**2/xy.shape[1]
s = directional.mean(axis=1)
frames = np.arange(1, len(s)+1)
assert np.isfinite(s).all() and (s >= 0).all()
previous = json.loads((a.analysis_dir/'summary.json').read_text())
np.testing.assert_allclose(s.mean(), previous['S_at_smallest_q'], rtol=1e-10)
with (a.analysis_dir/'smallest_q_by_frame.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['production_frame','timestep','q_sigma_inverse','S_qmin_shell','S_qx','S_qy'])
    w.writerows(zip(frames, steps, np.full(len(s), qmin), s, directional[:,0], directional[:,1]))
stats = dict(q_min=float(qmin), frames=len(s), mean=float(s.mean()),
    frame_standard_deviation=float(s.std(ddof=1)),minimum=float(s.min()),maximum=float(s.max()),
    minimum_frame=int(frames[s.argmin()]),maximum_frame=int(frames[s.argmax()]),
    block_means=[float(b.mean()) for b in np.array_split(s,5)],
    definition='Mean of S(qmin,0) and S(0,qmin) for each frame; not a radial bin.',
    caveat='Frames can be correlated. Frame standard deviation is not uncertainty in the mean.')
(a.analysis_dir/'smallest_q_summary.json').write_text(json.dumps(stats,indent=2))
plt.rcParams.update({'font.size':11,'axes.spines.right':False})
fig, ax = plt.subplots(figsize=(10,4.8))
for i, inds in enumerate(np.array_split(np.arange(len(s)),5)):
    if i%2 == 0:
        ax.axvspan(frames[inds[0]]-.5,frames[inds[-1]]+.5,color='#007c91',alpha=.045)
ax.plot(frames,s,'o-',lw=1.3,ms=3.2,color='#007c91',label='Individual frame (x/y shell average)')
ax.axhline(s.mean(),color='#424242',ls='--',lw=1.2,label=f'Overall mean = {s.mean():.5f}')
for i, inds in enumerate(np.array_split(np.arange(len(s)),5)):
    ax.hlines(s[inds].mean(),frames[inds[0]]-.5,frames[inds[-1]]+.5,
              color='#dd8525',lw=2.5,label='20-frame block mean' if i==0 else None)
ax.set(xlim=(.5,len(s)+.5),ylim=(0,float(s.max())*1.4),xlabel='Production frame',
       ylabel=r'$S_{\mathrm{core}}(q_{\min})$',
       title=rf'Smallest-wavevector fluctuations: $q_{{\min}}\sigma = {qmin:.5f}$')
ax.grid(alpha=.18)
ax.legend(loc='upper left',fontsize=9,framealpha=.95)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
fig.text(.5,.015,'400 cores · 100 frames · steps 2,610,000–3,600,000 · one frame every 10,000 steps',
         ha='center',fontsize=10)
fig.tight_layout(rect=[0,.045,1,1])
for ext in ['png','pdf']:
    fig.savefig(a.analysis_dir/f'smallest_q_by_frame.{ext}',dpi=220)
plt.close(fig)
print(json.dumps(stats,indent=2))
