import numpy as np, time
from pcm1d import *
pcm = pcm_properties('erythritol')
props = composite_properties(pcm, [0.0], 0.0)
Tw, Ti = pcm['Tm']+20, pcm['Tm']-10
alpha = pcm['k']/(pcm['rho']*pcm['cp'])
St_l = pcm['cp']*(Tw-pcm['Tm'])/pcm['L']; St_s = pcm['cp']*(pcm['Tm']-Ti)/pcm['L']
lam = neumann_lambda(St_l, St_s)
print('St_l',St_l,'St_s',St_s,'lambda',lam)
for N in [100,200,400]:
    g = make_grid('planar', 0.0, 0.10, N)
    src = dict(Tmean=Tw, A=np.zeros(1), P=np.ones(1), shape='constant', phase=np.zeros(1))
    t0=time.time()
    r = solve_pcm(props, g, 1e9, Ti, src, 3600*2, 2.0, record_every=600)
    # front: position where f crosses 0.5
    f = r['f'][0]; x = g['rc']
    i = np.argmax(f<0.5); xf = x[i-1] + (0.5-f[i-1])*(x[i]-x[i-1])/(f[i]-f[i-1])
    Ta, s = neumann_solution(x, 7200, Tw, Ti, pcm['Tm'], alpha, lam)
    err = np.max(np.abs(r['T'][0]-Ta))
    print(N, 'front num %.3f mm  analytic %.3f mm  rel err %.2f%%  maxTerr %.3f K  iters %.1f/%d  %.1fs'%(xf*1e3, s*1e3, 100*(xf-s)/s, err, r['mean_iter'], r['max_iter'], time.time()-t0))
    # energy balance
    print('   energy balance rel err', (r['E'][0,-1]-r['Qin'][0,-1])/r['Qin'][0,-1])
