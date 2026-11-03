import numpy as np, time
from pcm1d import *
from analysis import optimum_loading
from scipy.interpolate import PchipInterpolator
W=np.round(np.arange(0,0.2501,0.005),4); g=make_grid('radial',0.0125,0.025,80)
def run(pcm_name, ks, cases, dTm=1.0):
    pcm=dict(pcm_properties(pcm_name)); pcm['dTm']=dTm
    ww,cc=np.meshgrid(W,np.arange(len(cases)),indexing='ij'); ww=ww.ravel(); cc=cc.ravel(); C=np.array(cases)
    props=composite_properties(pcm,ww,K_SCENARIOS[ks])
    src=dict(Tmean=pcm['Tm']+20,shape_id=C[cc,0],A=C[cc,1],P=C[cc,2],phase=C[cc,3])
    r=solve_pcm(props,g,1000.0,pcm['Tm']-10,src,3600,2.0,record_every=60)
    return r,ww,cc
def curve(r,ww,cc,idx,tc):
    j=int(np.argmin(abs(r['t']-tc))); E=r['E'][:,j]
    Ew=[E[cc==q][np.argsort(ww[cc==q])] for q in idx]; return np.mean(Ew,0)
def metrics(Ec,Ef):
    wc,Ecs,_=optimum_loading(W,Ec); wf,Efs,_=optimum_loading(W,Ef)
    pen=100*(1-PchipInterpolator(W,Ef)(wc)/Efs); return wc,wf,pen
# ---- 1. mushy-zone width
print('MUSHY ZONE (erythritol & paraffin, low route, sine/square A=20 P=30, 4 phases)')
ph=[0,np.pi/2,np.pi,3*np.pi/2]
cases=[(0,0,3600,0)]+[(1,20,1800,p) for p in ph]+[(2,20,1800,p) for p in ph]
for pcm_name in ['erythritol','paraffin']:
    for dTm in [0.5,1.0,2.0]:
        t0=time.time(); r,ww,cc=run(pcm_name,'low',cases,dTm)
        out=[]
        for tc in [1800,3600]:
            Ec=curve(r,ww,cc,[0],tc); 
            for nm,idx in [('sine',[1,2,3,4]),('square',[5,6,7,8])]:
                wc,wf,pen=metrics(Ec,curve(r,ww,cc,idx,tc))
                out.append(f'tc{tc//60} {nm}: w*c={100*wc:.1f} w*f={100*wf:.1f} pen={pen:.2f}%')
        print(f'{pcm_name} dTm={dTm}: Econst60={curve(r,ww,cc,[0],3600).max()/1e3:.1f}  '+' | '.join(out), f'({time.time()-t0:.0f}s)', flush=True)
# ---- 2. 4 vs 8 phases (worst regime P=60 > tc=30, and tc=60)
print('PHASES (erythritol, low & mid, A=20, P=60)')
ph8=[k*np.pi/4 for k in range(8)]
cases=[(0,0,3600,0)]+[(1,20,3600,p) for p in ph8]+[(2,20,3600,p) for p in ph8]
for ks in ['low','mid']:
    r,ww,cc=run('erythritol',ks,cases)
    for tc in [1800,3600]:
        Ec=curve(r,ww,cc,[0],tc)
        for nm,base in [('sine',1),('square',9)]:
            i8=list(range(base,base+8)); i4=i8[::2]
            wc,wf8,p8=metrics(Ec,curve(r,ww,cc,i8,tc)); _,wf4,p4=metrics(Ec,curve(r,ww,cc,i4,tc))
            E8=curve(r,ww,cc,i8,tc).max(); E4=curve(r,ww,cc,i4,tc).max()
            print(f'{ks} tc{tc//60} {nm}: w* 4ph={100*wf4:.1f} 8ph={100*wf8:.1f} | penalty 4ph={p4:.2f}% 8ph={p8:.2f}% | E* diff {100*(E4/E8-1):+.2f}%', flush=True)
