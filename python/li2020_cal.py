import numpy as np
from pcm1d import *
pcm = pcm_properties('li2020_salt'); g = make_grid('radial', 0.0125, 0.025, 80)
hs = np.logspace(np.log10(20), np.log10(400), 40); B=len(hs)
props = composite_properties(pcm, np.zeros(B), 0.0)
src = dict(Tmean=300.0, A=np.zeros(B), P=np.full(B,3600.), shape_id=np.zeros(B), phase=np.zeros(B))
r = solve_pcm(props, g, hs, pcm['Tm']-10, src, 3*3600, 2.0, record_every=10)
tm = np.array([r['t'][np.argmax(r['fl'][b]>=0.999)]/60 if r['fl'][b].max()>=0.999 else np.nan for b in range(B)])
ok=~np.isnan(tm)
h_cal = np.interp(49.7, tm[ok][::-1], hs[ok][::-1]); print('calibrated h = %.1f W/m2K'%h_cal)
# now predict fluctuating cases with h_cal
cases=[('constant',0,3600),('P=1',100,60),('P=6',100,360),('P=20',100,1200),('P=60',100,3600),('A=50',50,3600),('A=150',150,3600)]
B=len(cases); props=composite_properties(pcm,np.zeros(B),0.0)
src=dict(Tmean=300.0,A=np.array([c[1] for c in cases],float),P=np.array([c[2] for c in cases],float),
         shape_id=np.array([0]+[1]*(B-1)),phase=np.zeros(B))
r=solve_pcm(props,g,h_cal,pcm['Tm']-10,src,3*3600,1.0,record_every=5)
li_t={'constant':49.7,'P=1':49.8,'P=6':49.3,'P=20':None,'P=60':37.3,'A=50':44.2,'A=150':32.3}
res=[]
for b,c in enumerate(cases):
    j=np.argmax(r['fl'][b]>=0.999); t=r['t'][j]/60; E=r['E'][b,j]/1e3
    res.append((c[0],t,E,li_t[c[0]]))
E0=res[0][2]; t0=res[0][1]
for c,t,E,lt in res:
    print('%-9s t_melt %.1f min (%+.1f%%)  Li %s   E_melt %.1f kJ/m (%+.1f%%)'%(c,t,100*(t/t0-1),lt,E,100*(E/E0-1)))
np.save('results/li2020_cal.npy',np.array([h_cal]+[x for r_ in res for x in r_],dtype=object),allow_pickle=True)
