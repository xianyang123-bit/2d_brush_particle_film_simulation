"""Compare core S(q) for single frames and consecutive 5/20-frame averages."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('analysis_dir', type=Path)
a = p.parse_args()
d = np.load(a.analysis_dir / 'structure_factor_frames.npz')
q, s, steps = d['q'], d['S_by_frame'], d['steps']
assert s.shape == (100, len(q)) and np.isfinite(s).all()
mean = s.mean(axis=0)
windows = {n: s.reshape(-1,n,len(q)).mean(axis=1) for n in (1,5,20)}
for n, curves in windows.items():
    np.testing.assert_allclose(curves.mean(axis=0), mean, atol=1e-12)
    with (a.analysis_dir / f'Sq_{n:02d}_frame_windows.csv').open('w',newline='') as f:
        w = csv.writer(f)
        w.writerow(['window','first_step','last_step','q_sigma_inverse','S_q'])
        for i, curve in enumerate(curves):
            w.writerows((i+1,int(steps[i*n]),int(steps[(i+1)*n-1]),float(x),float(y))
                        for x,y in zip(q,curve))

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig, axes = plt.subplots(2,3,figsize=(14,7.5),sharex='row',sharey='row',layout='constrained')
cmap = plt.get_cmap('viridis')
norm = Normalize(float(steps[0]/1e6),float(steps[-1]/1e6))
low = q <= .6
low_ymax = 1.08 * max(s[:,low].max(), max(np.interp(.6,q,curve) for curve in s))
for col,(n,curves) in enumerate(windows.items()):
    centers = steps.reshape(-1,n).mean(axis=1)/1e6
    for row in range(2):
        ax = axes[row,col]
        for curve,t in zip(curves,centers):
            ax.plot(q,curve,color=cmap(norm(t)),alpha={1:.24,5:.6,20:.9}[n],
                    lw={1:.65,5:.9,20:1.2}[n])
        ax.plot(q,mean,color='#171717',lw=1.5,ls='--',label='100-frame mean',zorder=5)
        ax.grid(alpha=.16)
        ax.set_xlabel(r'$q\sigma$')
        if col == 0:
            ax.set_ylabel(r'$S_{\mathrm{core}}(q)$')
        if row == 0:
            ax.axhline(1,color='gray',lw=.75,ls=':')
            ax.set(xlim=(0,6),ylim=(0,1.05*s.max()))
            ax.set_title(f'{n} frame'+('' if n==1 else 's')+f' per curve\n{len(curves)} separate curves')
        else:
            ax.set(xlim=(0,.6),ylim=(0,low_ymax))
axes[0,0].legend(loc='upper right',fontsize=8)
axes[1,0].text(.03,.93,'Small-q detail',transform=axes[1,0].transAxes,va='top',fontsize=10)
cb = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),pad=.015,shrink=.85)
cb.set_label('Window center timestep (millions)')
fig.suptitle('Core structure factor: effect of averaging consecutive frames\n'
             '400 cores · 100 production frames · identical q bins and axis scales across columns',fontsize=13)
fig.supxlabel('Nonoverlapping windows; frames are 10,000 steps apart. Temporal correlations may persist.',fontsize=10)
for ext in ('png','pdf'):
    fig.savefig(a.analysis_dir/f'frame_averaging_comparison.{ext}',dpi=220)
plt.close(fig)

stds = {n:v.std(axis=0,ddof=1) for n,v in windows.items()}
with (a.analysis_dir/'frame_averaging_variation.csv').open('w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['q_sigma_inverse','mean_S_q','std_1_frame','std_5_frames','std_20_frames'])
    w.writerows(zip(q,mean,stds[1],stds[5],stds[20]))
peak=int(np.argmax(mean))
summary={'frames':len(steps),'first_step':int(steps[0]),'last_step':int(steps[-1]),
         'frame_stride_steps':int(steps[1]-steps[0]),
         'peak_q':float(q[peak]),'peak_mean':float(mean[peak]),
         'window_comparison':[{'frames_per_window':n,'windows':len(v),
             'peak_standard_deviation':float(stds[n][peak]),
             'peak_min':float(v[:,peak].min()),'peak_max':float(v[:,peak].max())}
             for n,v in windows.items()],
         'interpretation':'Standard deviation across nonoverlapping window curves, not uncertainty in the overall mean. Windows can be correlated. Only five 20-frame windows are available.'}
(a.analysis_dir/'frame_averaging_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
