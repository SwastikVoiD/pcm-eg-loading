import numpy as np, time
from pcm1d import *
from analysis import optimum_loading
from scipy.interpolate import PchipInterpolator
W=np.round(np.arange(0,0.2501,0.005),4); g=make_grid('radial',0.0125,0.025,80)
ph16=[k*np.pi/8 for k in range(16)]
cases=[(0,0,3600,0)]
for sid in [1,2]:
    for A in [10,20]:
        cases += [(sid,A,3600,p) for p in ph16]
C=np.array(cases); rows=[]
for pcm_name in ['erythritol','paraffin']:
  for ks in ['low','mid','high']:
    pcm=pcm_properties(pcm_name); t0=time.time()
    ww,cc=np.meshgrid(W,np.arange(len(cases)),indexing='ij'); ww=ww.ravel(); cc=cc.ravel()
    props=composite_properties(pcm,ww,K_SCENARIOS[ks])
    src=dict(Tmean=pcm['Tm']+20,shape_id=C[cc,0],A=C[cc,1],P=C[cc,2],phase=C[cc,3])
    r=solve_pcm(props,g,1000.0,pcm['Tm']-10,src,1800,2.0,record_every=60)
    E=r['E'][:,-1]
    cur=lambda idx: np.mean([E[cc==q][np.argsort(ww[cc==q])] for q in idx],0)
    Ec=cur([0]); wc,Ecs,_=optimum_loading(W,Ec)
    for gi,(sid,A) in enumerate([(1,10),(1,20),(2,10),(2,20)]):
        i16=list(range(1+16*gi,17+16*gi))
        for nph,idx in [(4,i16[::4]),(8,i16[::2]),(16,i16)]:
            Ef=cur(idx); wf,Efs,_=optimum_loading(W,Ef)
            pen=100*(1-PchipInterpolator(W,Ef)(wc)/Efs); dE=100*(Efs/Ecs-1)
            rows.append((pcm_name,ks,['','sine','square'][sid],A,nph,wc,wf,wf-wc,pen,dE))
    print(pcm_name,ks,f'{time.time()-t0:.0f}s',flush=True)
import pandas as pd
D=pd.DataFrame(rows,columns=['pcm','k','shape','A','nph','w_const','w_fluct','dw','penalty','dE'])
D.to_csv('results/phase_convergence_Pgt_tc.csv',index=False)
pd.set_option('display.width',200)
print(D.pivot_table(index=['pcm','k','shape','A'],columns='nph',values=['penalty','dw','dE']).round(3).to_string())
for n in [4,8,16]:
    s=D[D.nph==n]; print(n,'phases: max penalty %.2f%%, max |dw| %.3f, min dE %.1f%%'%(s.penalty.max(),s.dw.abs().max(),s.dE.min()))
