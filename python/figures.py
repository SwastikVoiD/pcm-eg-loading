import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pcm1d import *
from analysis import summarize, optimum_loading
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.3, 'figure.dpi': 200,
                     'font.family': 'DejaVu Sans'})
C3 = ['#1b6ca8', '#d95f02', '#2a9d4b']
OUT = 'figures'

def fig1_k_models():
    w = np.linspace(0.002, 0.25, 200)
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for (nm, b), c in zip(K_SCENARIOS.items(), C3):
        ax.plot(100*w, b*w, color=c, lw=1.6, label=f'{nm}: β = {b:.1f} W m$^{{-1}}$K$^{{-1}}$')
    e = pcm_properties('erythritol'); phi = (w/2200)/((w/2200)+(1-w)/e['rho'])
    ax.plot(100*w, k_maxwell(0.72, 150.0, phi) - 0.72, 'k--', lw=1.2, label='Maxwell, dispersed spheres')
    ax.plot([2,4,7,10], np.array([.40,.52,.68,.82])-.22, 'o', color=C3[0], label='Sari & Karaipekli 2007 (paraffin)')
    ax.plot([20], [6.003-.348], 's', color=C3[1], label='Li et al. 2024 (paraffin)')
    ax.plot([10], [12.51-0.72], '^', color=C3[2], label='Yuan et al. 2019 (erythritol)')
    ax.set_yscale('log'); ax.set_ylim(0.01, 50); ax.set_xlim(0, 25)
    ax.set_xlabel('EG mass fraction w (%)'); ax.set_ylabel('Conductivity gain k − k$_{pcm}$ (W m$^{-1}$ K$^{-1}$)')
    ax.legend(fontsize=6.5, loc='lower right'); fig.tight_layout(); fig.savefig(f'{OUT}/fig1_k_models.png'); plt.close(fig)

def fig2_validation():
    pcm = pcm_properties('erythritol'); props = composite_properties(pcm, [0.0], 0.0)
    Tw, Ti = pcm['Tm']+20, pcm['Tm']-10; alpha = pcm['k']/(pcm['rho']*pcm['cp'])
    lam = neumann_lambda(pcm['cp']*20/pcm['L'], pcm['cp']*10/pcm['L'])
    g = make_grid('planar', 0, 0.3, 600)
    src = dict(Tmean=Tw, A=np.zeros(1), P=np.ones(1), shape_id=np.zeros(1), phase=np.zeros(1))
    out = {}
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for tt, c in zip([1800, 7200], C3):
        r = solve_pcm(props, g, 1e9, Ti, src, tt, 2.0, record_every=600)
        Ta, s = neumann_solution(g['rc'], tt, Tw, Ti, pcm['Tm'], alpha, lam)
        m = g['rc'] < 0.06
        ax.plot(1e3*g['rc'][m], Ta[m], '-', color=c, lw=1.4, label=f'Neumann, t = {tt//60} min')
        ax.plot(1e3*g['rc'][m][::6], r['T'][0][m][::6], 'o', ms=3.2, mfc='none', color=c, label=f'Model, t = {tt//60} min')
        f = r['f'][0]; x = g['rc']; i = np.argmax(f < 0.5)
        xf = x[i-1] + (0.5-f[i-1])*(x[i]-x[i-1])/(f[i]-f[i-1])
        e = r['T'][0] - Ta
        out[tt] = dict(front_num=1e3*xf, front_exact=1e3*s, err=100*(xf-s)/s, maxT=np.max(np.abs(e[x < 0.08])),
                       bal=(r['E'][0,-1]-r['Qin'][0,-1])/r['Qin'][0,-1])
    ax.axhline(pcm['Tm'], color='k', ls=':', lw=0.8)
    ax.set_xlabel('Distance from heated wall (mm)'); ax.set_ylabel('Temperature (°C)')
    ax.legend(fontsize=6.5); fig.tight_layout(); fig.savefig(f'{OUT}/fig2_validation.png'); plt.close(fig)
    return out, lam

def fig3_transient():
    pcm = pcm_properties('erythritol'); props = composite_properties(pcm, [0.05]*3, K_SCENARIOS['low'])
    src = dict(Tmean=pcm['Tm']+20, A=np.array([0, 20., 20.]), P=np.full(3, 3600.), shape_id=np.array([0, 1, 2]), phase=np.zeros(3))
    g = make_grid('radial', 0.0125, 0.025, 80)
    r = solve_pcm(props, g, 1000.0, pcm['Tm']-10, src, 7200, 2.0, record_every=30)
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(4.6, 5.6), sharex=True)
    for b, (nm, c) in enumerate(zip(['constant', 'sine', 'square'], C3)):
        a1.plot(r['t']/60, r['Thtf'][b], color=c, lw=1.3, label=nm)
        a2.plot(r['t']/60, r['fl'][b], color=c, lw=1.3)
        a3.plot(r['t']/60, r['E'][b]/1e3, color=c, lw=1.3)
    a1.axhline(pcm['Tm'], color='k', ls=':', lw=0.8); a1.set_ylabel('HTF temp. (°C)'); a1.legend(fontsize=7, ncol=3, loc='lower left')
    a2.set_ylabel('Melt fraction (–)'); a3.set_ylabel('Stored energy (kJ m$^{-1}$)'); a3.set_xlabel('Time (min)')
    for a in (a1, a2, a3): a.axvline(30, color='grey', lw=0.6, ls='--'); a.axvline(60, color='grey', lw=0.6, ls='--')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig3_transient_example.png'); plt.close(fig)
    return r

def phase_avg_curve(d, t_c, key):
    j = int(np.argmin(np.abs(d['t'] - t_c))); E = d['E'][:, j]
    cidx = [i for i in range(len(d['cases'])) if np.allclose(d['cases'][i, :3], key)]
    w = np.unique(d['w']); Ew = []
    for q in cidx:
        sel = d['case'] == q; o = np.argsort(d['w'][sel]); Ew.append(E[sel][o])
    return w, np.mean(Ew, 0)

def fig4_E_vs_w(pcm_name='erythritol', ks='low', h=1000):
    """Each panel uses a source period equal to its window (P = t_c), the regime
    where four start phases are converged."""
    d = dict(np.load(f'results/sweep_{pcm_name}_{ks}_h{h}.npz'))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1)); out = {}
    for ax, tc in zip(axs, [1800, 3600]):
        P = tc
        for key, nm, c in zip([(0,0,3600),(1,20,P),(2,20,P)], ['constant', f'sine, A=20 K, P={P//60} min', f'square, A=20 K, P={P//60} min'], C3):
            w, Eb = phase_avg_curve(d, tc, key)
            wo, Eo, _ = optimum_loading(w, Eb); out[(tc, key[0])] = (wo, Eo)
            ax.plot(100*w, Eb/1e3, color=c, lw=1.5, label=nm); ax.plot(100*wo, Eo/1e3, 'o', color=c, ms=5)
        ax.axvline(10, color='k', ls=':', lw=0.8)
        ax.set_title(f'Charging window {tc//60} min', fontsize=9); ax.set_xlabel('EG mass fraction w (%)')
        ax.legend(fontsize=6.5, loc='lower center')
    axs[0].set_ylabel('Stored energy (kJ m$^{-1}$)')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig4_E_vs_w_{pcm_name}_{ks}.png'); plt.close(fig)
    return out

def fig5_optimum_map(ALL, h=1000):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
    for ax, pcm_name in zip(axs, ['erythritol', 'paraffin']):
        for i, (ks, c) in enumerate(zip(['low', 'mid', 'high'], C3)):
            T = ALL[(ALL.pcm == pcm_name) & (ALL.k_scenario == ks) & (ALL.h == h)]
            tcs = sorted(T.t_c_min.unique()); wc = []; lo = []; hi = []
            for tc in tcs:
                r = T[T.t_c_min == tc]
                wc.append(r[r['shape'] == 'constant'].w_opt.iloc[0]); lo.append(r.w_opt.min()); hi.append(r.w_opt.max())
            wc, lo, hi = map(np.array, (wc, lo, hi)); x = np.array(tcs) + (i-1)*2.5
            ax.errorbar(x, 100*wc, yerr=[100*(wc-lo), 100*(hi-wc)], fmt='-o', color=c, ms=4, capsize=3, lw=1.3, label=f'{ks}-k route')
        ax.axhline(10, color='k', ls=':', lw=0.8); ax.set_title(pcm_name, fontsize=9)
        ax.set_xlabel('Charging window (min)'); ax.set_xticks([30, 60, 120])
    axs[0].set_ylabel('Optimum EG fraction w* (%)'); axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig5_optimum_map.png'); plt.close(fig)

def fig6_penalty(ALL):
    T = ALL[(ALL['shape'] != 'constant') & (~ALL.bypass) & (ALL.h == 1000)].copy()
    T['ratio'] = T.P_min / T.t_c_min
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    rng = np.random.default_rng(0)
    for (pcm_name, mk) in [('erythritol', 'o'), ('paraffin', 's')]:
        for ks, c in zip(['low', 'mid', 'high'], C3):
            r = T[(T.pcm == pcm_name) & (T.k_scenario == ks)]
            x = r.ratio * np.exp(rng.uniform(-0.08, 0.08, len(r)))
            lab = f'{pcm_name}, {ks}'
            a1.scatter(x, r.penalty_pct, s=10, marker=mk, color=c, alpha=0.6, label=lab)
            a2.scatter(x, r.dE_pct, s=10, marker=mk, color=c, alpha=0.6)
    for a in (a1, a2):
        a.set_xscale('log'); a.axvline(1, color='k', ls=':', lw=0.8); a.set_xlabel('Source period / charging window, P/t$_c$')
        a.set_xticks([1/12, 1/4, 1/2, 1, 2]); a.set_xticklabels(['1/12', '1/4', '1/2', '1', '2'])
    a1.set_ylabel('Design penalty of steady-source w* (%)'); a2.set_ylabel('Change in achievable energy (%)')
    a1.legend(fontsize=6, ncol=2, loc='upper left')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig6_penalty_vs_ratio.png'); plt.close(fig)

def fig7_bypass(ALL, tc=60):
    T = ALL[(ALL.h == 1000) & (ALL.t_c_min == tc) & (ALL['shape'] != 'constant')]
    groups = [('erythritol', 'low'), ('erythritol', 'mid'), ('erythritol', 'high'), ('paraffin', 'mid')]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    rng = np.random.default_rng(1)
    for i, (p, k) in enumerate(groups):
        for bp, c, dx, lab in [(False, C3[0], -0.15, 'HTF always on'), (True, C3[1], 0.15, 'bypass when cold')]:
            r = T[(T.pcm == p) & (T.k_scenario == k) & (T.bypass == bp)]
            x = i + dx + rng.uniform(-0.06, 0.06, len(r))
            a1.scatter(x, r.dE_pct, s=10, color=c, alpha=0.7, label=lab if i == 0 else None)
            a2.scatter(x, r.penalty_pct, s=10, color=c, alpha=0.7)
    for a in (a1, a2):
        a.set_xticks(range(4)); a.set_xticklabels([f'{p[:4]}.\n{k}' for p, k in groups])
    a1.axhline(0, color='k', lw=0.8, ls=':')
    a1.set_ylabel('Change in achievable energy (%)'); a2.set_ylabel('Design penalty of steady-source w* (%)')
    a1.legend(fontsize=7); fig.suptitle(f'Charging window {tc} min, h = 1000 W m$^{{-2}}$ K$^{{-1}}$', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig7_bypass.png'); plt.close(fig)

def fig8_convection_bound(pcm_name='erythritol', ks='low'):
    d0 = dict(np.load(f'results/sweep_{pcm_name}_{ks}_h1000.npz'))
    d1 = dict(np.load(f'results/sweep_{pcm_name}_{ks}_h1000_conv.npz'))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, tc in zip(axs, [1800, 3600]):
        for d, lab, ls in [(d0, 'conduction only (lower bound)', '-'), (d1, 'with convection (upper bound)', '--')]:
            for key, nm, c in [((0, 0, 3600), 'constant', C3[0]), ((2, 20, 1800), 'square, A=20 K, P=30 min', C3[2])]:
                w, Eb = phase_avg_curve(d, tc, key)
                wo, Eo, _ = optimum_loading(w, Eb)
                ax.plot(100*w, Eb/1e3, ls, color=c, lw=1.4, label=f'{nm}, {lab}')
                ax.plot(100*wo, Eo/1e3, 'o', color=c, ms=4, mfc='white' if ls == '--' else c)
        ax.axvline(10, color='k', ls=':', lw=0.8)
        ax.set_title(f'Charging window {tc//60} min', fontsize=9); ax.set_xlabel('EG mass fraction w (%)')
    axs[0].set_ylabel('Stored energy (kJ m$^{-1}$)'); axs[0].legend(fontsize=5.8, loc='lower right')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig8_convection_bound_{pcm_name}_{ks}.png'); plt.close(fig)
