import numpy as np, pandas as pd, glob, re
from analysis import summarize
rows=[]
for fn in sorted(glob.glob('results/sweep_*.npz')):
    m=re.match(r'results/sweep_(\w+?)_(low|mid|high)_h(\d+)(_bypass)?(_conv)?\.npz',fn)
    pcm,ks,h,bp,cv=m.group(1),m.group(2),int(m.group(3)),bool(m.group(4)),bool(m.group(5))
    d=dict(np.load(fn))
    for tc in [1800,3600,7200]:
        S,_,_=summarize(d,tc); S.insert(0,'convection',cv); S.insert(0,'bypass',bp); S.insert(0,'t_c_min',tc//60); S.insert(0,'h',h); S.insert(0,'k_scenario',ks); S.insert(0,'pcm',pcm)
        rows.append(S)
ALL=pd.concat(rows,ignore_index=True); ALL.to_csv('results/summary_all.csv',index=False)
print(len(ALL))
