import numpy as np, time
from pcm1d import *
pcm = pcm_properties('erythritol')
w = np.array([0.0, 0.05, 0.10, 0.20])
props = composite_properties(pcm, w, K_SCENARIOS['mid'])
B=len(w); Tm=pcm['Tm']
src=dict(Tmean=Tm+20,A=np.full(B,20.0),P=np.full(B,1800.0),shape='sine',phase=np.zeros(B))
res={}
for N,dt in [(40,2.0),(80,2.0),(160,2.0),(80,1.0),(80,5.0)]:
    g=make_grid('radial',0.0125,0.025,N); t0=time.time()
    r=solve_pcm(props,g,1000.0,Tm-10,src,3600,dt,record_every=60)
    res[(N,dt)]=r['E'][:,-1]
    print(N,dt,'E(60min) kJ/m:',np.round(r['E'][:,-1]/1e3,3),'fl:',np.round(r['fl'][:,-1],4),'%.1fs'%(time.time()-t0), 'iters %.2f max %d'%(r['mean_iter'],r['max_iter']))
ref=res[(160,2.0)]
for key in res: print(key,'max rel dev vs N=160 dt=2: %.4f%%'%(100*np.max(np.abs(res[key]-ref)/ref)))
ref=res[(80,1.0)]
for key in [(80,2.0),(80,5.0)]: print(key,'max rel dev vs dt=1: %.4f%%'%(100*np.max(np.abs(res[key]-ref)/ref)))
