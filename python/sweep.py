import numpy as np, time, itertools, sys
from pcm1d import *

W = np.round(np.arange(0, 0.2501, 0.005), 4)
PERIODS = [600.0, 1800.0, 3600.0]; AMPS = [10.0, 20.0]
PHASES = [0.0, np.pi/2, np.pi, 3*np.pi/2]
DT_SUPER, DT_SUB = 20.0, 10.0
R_IN, R_OUT, NCELL, DT, T_END = 0.0125, 0.025, 80, 2.0, 7200.0

def build_cases():
    cases = [(0, 0.0, 3600.0, 0.0)]   # constant
    for sid, P, A, ph in itertools.product([1, 2], PERIODS, AMPS, PHASES):
        cases.append((sid, A, P, ph))
    return cases

def run(pcm_name, kscen, h, bypass=False, conv=False):
    pcm = pcm_properties(pcm_name); cases = build_cases()
    ww, cc = np.meshgrid(W, np.arange(len(cases)), indexing='ij')
    ww = ww.ravel(); cc = cc.ravel()
    props = composite_properties(pcm, ww, K_SCENARIOS[kscen])
    C = np.array(cases)
    src = dict(Tmean=pcm['Tm'] + DT_SUPER, shape_id=C[cc, 0], A=C[cc, 1], P=C[cc, 2], phase=C[cc, 3])
    g = make_grid('radial', R_IN, R_OUT, NCELL)
    t0 = time.time()
    cv = convection_params(pcm_name, ww) if conv else None
    r = solve_pcm(props, g, h, pcm['Tm'] - DT_SUB, src, T_END, DT, record_every=60.0, bypass=bypass, conv=cv)
    np.savez_compressed(f'results/sweep_{pcm_name}_{kscen}_h{int(h)}' + ('_bypass' if bypass else '') + ('_conv' if conv else '') + '.npz', t=r['t'], E=r['E'], fl=r['fl'], Qin=r['Qin'],
                        w=ww, case=cc, cases=C, k=props['k'], rho=props['rho'], L=props['L'], cp=props['cp'],
                        mean_iter=r['mean_iter'], max_iter=r['max_iter'])
    bal = np.max(np.abs(r['E'][:, -1] - r['Qin'][:, -1]) / np.abs(r['Qin'][:, -1]))
    print(pcm_name, kscen, h, 'bypass' if bypass else '', 'conv' if conv else '', 'B=%d  %.0fs  iters %.2f/%d  energy-balance max rel err %.2e' % (len(ww), time.time()-t0, r['mean_iter'], r['max_iter'], bal), flush=True)

if __name__ == '__main__':
    import os; os.makedirs('results', exist_ok=True)
    jobs = [(p, k, h) for p in ['erythritol', 'paraffin'] for k in ['low', 'mid', 'high'] for h in [1000.0, 200.0, 5000.0]]
    if len(sys.argv) > 1: jobs = jobs[:int(sys.argv[1])]
    for j in jobs: run(*j)
