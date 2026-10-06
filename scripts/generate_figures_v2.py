#!/usr/bin/env python3
"""
generate_figures_v2.py -- Figures for Paper III, version 2.
Reads data/singular_series_exact.csv, data/singular_series_convergence_v2.csv,
data/quadruplets_742.csv.  Writes figures/v2/v2_fig1.png ... v2_fig5.png (+pdf).

Usage:  python generate_figures_v2.py
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, optimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'figures', 'v2')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.size': 11, 'axes.labelsize': 12,
                     'axes.titlesize': 12.5, 'figure.dpi': 150, 'savefig.dpi': 300,
                     'savefig.bbox': 'tight', 'mathtext.fontset': 'cm'})
LN47 = math.log(47)

S = {}
with open(os.path.join(DATA, 'singular_series_exact.csv')) as f:
    for r in csv.DictReader(f):
        S[int(r['k'])] = float(r['S_k'])
quads = np.array(sorted(int(r['starting_n']) for r in csv.DictReader(open(os.path.join(DATA, 'quadruplets_742.csv')))))
quints = np.array([35676017721, 64482563907, 73292417435, 116255850744,
                   147743683226, 159430471996, 182501065420])


def I(k, N):
    return integrate.quad(lambda t: 1.0 / (46 * math.log(t) + LN47) ** k, 2, N, limit=500)[0]


def Q_mod(n, p):
    return (pow(n, 47, p) - pow(n - 1, 47, p)) % p


def sieve(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png')); fig.savefig(os.path.join(OUT, name + '.pdf'))
    plt.close(fig); print('  saved', name)


# ------------------------------------------------------------ Figure 1
def figure1():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), height_ratios=[1.6, 1.6],
                                   gridspec_kw={'hspace': 0.35})
    ax1.scatter(quads / 1e9, np.ones_like(quads), c='#9a9a9a', s=14, alpha=0.5, marker='|')
    for i, q in enumerate(quints):
        ax1.scatter(q / 1e9, 1.0, c='gold', s=260, marker='*', edgecolors='darkorange', zorder=5)
        ax1.annotate(f'Q{i+1}', (q / 1e9, 1.0), xytext=(q / 1e9, 1.22), fontsize=8,
                     ha='center', color='#8B4513', fontweight='bold')
    ax1.set_xlim(-5, 210); ax1.set_ylim(0.6, 1.5); ax1.set_yticks([])
    ax1.set_xlabel(r'$n$ (billions)')
    ax1.set_title('(a) The 7 quintuplets (stars) within the field of 742 quadruplets (ticks)', fontweight='bold')

    Nq = np.sort(quints); C = np.arange(1, 8)
    ax2.step(np.concatenate([[0], Nq / 1e9, [200]]), np.concatenate([[0], C, [7]]),
             color='darkorange', lw=2.5, where='post', label='Observed $C_5(N)$')
    Nf = np.linspace(1e9, 2.05e11, 300)
    bh = [S[5] * I(5, n) for n in Nf]
    ax2.plot(Nf / 1e9, bh, '--', color='#27ae60', lw=1.6,
             label=r'Bateman–Horn, $\mathfrak{S}_5 = %s$ (exact, no calibration)' % f"{S[5]:,.0f}")
    ax2.fill_between(Nf / 1e9, np.array(bh) - np.sqrt(bh), np.array(bh) + np.sqrt(bh),
                     color='#27ae60', alpha=0.12, label=r'$\pm 1\sigma$ (Poisson)')
    ax2.set_xlabel(r'$N$ (billions)'); ax2.set_ylabel(r'$C_5(N)$'); ax2.set_xlim(0, 210); ax2.set_ylim(0, 10)
    ax2.legend(fontsize=9.5, loc='upper left'); ax2.grid(True, alpha=0.2)
    ax2.set_title('(b) Cumulative quintuplet count', fontweight='bold')
    save(fig, 'v2_fig1')


# ------------------------------------------------------------ Figure 2
def figure2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    labels = ['k=1\n(2e9)', 'k=2\n(2e9)', 'k=3\n(2e9)', 'k=4\n(2e9)', 'k=4\n(2e11)', 'k=5\n(2e11)', 'k=6\n(2e11)']
    # observed POSITION counts (maximal runs converted: a run of length m gives m-k+1 positions)
    obs = [18121562 + 2 * 173351 + 3 * 1749 + 4 * 15, 173351 + 2 * 1749 + 3 * 15, 1749 + 2 * 15, 15,
           735 + 2 * 7, 7, 0]
    pred = [S[1] * I(1, 2e9), S[2] * I(2, 2e9), S[3] * I(3, 2e9), S[4] * I(4, 2e9),
            S[4] * I(4, 2e11), S[5] * I(5, 2e11), S[6] * I(6, 2e11)]
    x = np.arange(7); w = 0.38
    ax1.bar(x - w / 2, [math.log10(max(o, 0.03)) for o in obs], w, color='#e74c3c', alpha=0.85, label='Observed (positions)')
    ax1.bar(x + w / 2, [math.log10(p) for p in pred], w, color='#3498db', alpha=0.85, label='Bateman–Horn (exact $\\mathfrak{S}_k$)')
    for i, (o, p) in enumerate(zip(obs, pred)):
        lo = f'{o:,}' if o else '0'
        lp = f'{p:,.0f}' if p >= 10 else f'{p:.2f}'
        if o >= 10:
            ax1.text(i - w / 2, math.log10(o) - 0.12, lo, ha='center', va='top', rotation=90, fontsize=7, color='white', fontweight='bold')
            ax1.text(i + w / 2, math.log10(p) - 0.12, lp, ha='center', va='top', rotation=90, fontsize=7, color='white')
        else:
            ax1.text(i - w / 2, max(math.log10(max(o, 0.03)), 0) + 0.12, lo, ha='center', fontsize=7.5, fontweight='bold')
            ax1.text(i + w / 2, max(math.log10(p), 0) + 0.12, lp, ha='center', fontsize=7.5)
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=8.5)
    ax1.set_ylabel(r'$\log_{10}$ count'); ax1.set_ylim(-1.8, 8.2); ax1.grid(True, alpha=0.15, axis='y')
    ax1.set_title('(a) Constellation hierarchy: observed vs predicted', fontweight='bold'); ax1.legend(fontsize=9)

    n = np.linspace(1e10, 2e11, 300)
    p5 = (S[5] / S[4]) / (46 * np.log(n) + LN47); p6 = (S[6] / S[5]) / (46 * np.log(n) + LN47)
    ax2.plot(n / 1e9, p5 * 100, color='#9b59b6', lw=2, label=r'$P(\mathrm{quint}\mid\mathrm{quad})$')
    ax2.plot(n / 1e9, p6 * 100, color='#e74c3c', lw=2, label=r'$P(\mathrm{sext}\mid\mathrm{quint})$')
    for q in quints:
        ax2.plot(q / 1e9, 100 * (S[5] / S[4]) / (46 * math.log(q) + LN47), '*', color='gold', ms=12,
                 markeredgecolor='darkorange', zorder=5)
    ax2.set_xlabel(r'$n$ (billions)'); ax2.set_ylabel('Conditional probability (%)'); ax2.set_xlim(0, 205)
    ax2.set_title('(b) Extension probabilities', fontweight='bold'); ax2.legend(fontsize=10); ax2.grid(True, alpha=0.2)
    fig.tight_layout(); save(fig, 'v2_fig2')


# ------------------------------------------------------------ Figure 3
def figure3():
    fig, ax = plt.subplots(figsize=(10, 6))
    N = np.logspace(10, 14, 300)
    E5 = [S[5] * I(5, x) for x in N]; E6 = [S[6] * I(6, x) for x in N]
    ax.loglog(N, E6, color='#e74c3c', lw=2.5, label=r'$E[C_6(N)]$, $\mathfrak{S}_6 = %s$' % f"{S[6]:,.0f}")
    ax.loglog(N, E5, color='#9b59b6', lw=2, label=r'$E[C_5(N)]$, $\mathfrak{S}_5 = %s$' % f"{S[5]:,.0f}")
    ax.axvline(2e11, color='#3498db', ls='--', alpha=0.7); ax.text(2.3e11, 0.005, 'Current limit\n$N=2\\times10^{11}$', fontsize=9, color='#3498db')
    ax.axhline(1, color='grey', ls=':', alpha=0.6); ax.text(1.2e10, 1.3, r'$E[C_6]=1$', fontsize=10, color='grey')
    ax.plot(2e11, 7, 'o', color='#9b59b6', ms=10, markeredgecolor='k', zorder=5, label='Observed: 7 quintuplets')
    ax.plot(2e11, S[6] * I(6, 2e11), 's', color='#e74c3c', ms=10, markeredgecolor='k', zorder=5, label='Observed: 0 sextuplets (plotted at $E$)')
    Ns = 10 ** optimize.brentq(lambda lg: S[6] * I(6, 10 ** lg) - 1, 12, 15)
    ax.plot(Ns, 1, 'D', color='#e74c3c', ms=10, markeredgecolor='k', zorder=5)
    ax.annotate(f'$N^* \\approx {Ns/1e13:.2f}\\times10^{{13}}$', (Ns, 1), xytext=(Ns * 2.2, 0.3), fontsize=10,
                color='#c0392b', fontweight='bold', arrowprops=dict(arrowstyle='->', color='#c0392b'))
    ax.set_xlabel(r'Search bound $N$'); ax.set_ylabel('Expected count'); ax.set_xlim(5e9, 1.5e14); ax.set_ylim(3e-3, 500)
    ax.set_title('Predicted quintuplet and sextuplet counts (exact singular series)', fontweight='bold')
    ax.legend(fontsize=9.5, loc='upper left'); ax.grid(True, alpha=0.2, which='both')
    save(fig, 'v2_fig3')


# ------------------------------------------------------------ Figure 4
def figure4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ps = sieve(50); f = {4: [], 5: [], 6: []}
    for p in ps:
        for k in f:
            om = p - sum(1 for r in range(p) if all(Q_mod(r + i, p) for i in range(k)))
            f[k].append((1 - om / p) / (1 - 1 / p) ** k)
    x = np.arange(len(ps)); w = 0.25
    ax1.bar(x - w, f[4], w, color='#e74c3c', alpha=.8, label='$k=4$'); ax1.bar(x, f[5], w, color='#9b59b6', alpha=.8, label='$k=5$')
    ax1.bar(x + w, f[6], w, color='#3498db', alpha=.8, label='$k=6$')
    ax1.set_xticks(x); ax1.set_xticklabels(ps, fontsize=8); ax1.set_xlabel('Prime $p$'); ax1.set_ylabel(r'local factor $\sigma_k(p)$')
    ax1.set_title(r'(a) Local factors, $p<50$ (all non-resonant, $\omega_k=0$)', fontweight='bold'); ax1.legend(fontsize=9); ax1.grid(True, alpha=.2, axis='y')
    with open(os.path.join(DATA, 'singular_series_convergence_v2.csv')) as fh:
        rows = [r for r in csv.DictReader(fh)]
    allp = sieve(2000); cum = {4: [], 5: [], 6: []}; run = {4: 1.0, 5: 1.0, 6: 1.0}
    for p in allp:
        for k in cum:
            om = p - sum(1 for r in range(p) if all(Q_mod(r + i, p) for i in range(k))) if p % 47 == 1 else 0
            run[k] *= (1 - om / p) / (1 - 1 / p) ** k; cum[k].append(run[k])
    for k, c in [(4, '#e74c3c'), (5, '#9b59b6'), (6, '#3498db')]:
        ax2.semilogy(allp, cum[k], color=c, lw=2, label=r'$\mathfrak{S}_%d(B)$' % k)
        ax2.axhline(S[k], color=c, ls='--', lw=1, alpha=0.7)
    ax2.text(1250, S[4] * 1.15, 'dashed: exact limits', fontsize=8, color='#555')
    for p in (283, 659, 941, 1129, 1223, 1693, 1787):
        ax2.axvline(p, color='#333', ls=':', alpha=0.35)
    ax2.text(290, 2.5e6, '$p=283$\n(first cliff)', fontsize=8)
    ax2.set_xlabel('Sieve bound $B$'); ax2.set_ylabel(r'partial product $\mathfrak{S}_k(B)$'); ax2.set_ylim(1e3, 5e6)
    ax2.set_title('(b) Euler product with the periodic obstruction', fontweight='bold'); ax2.legend(fontsize=10, loc='lower right'); ax2.grid(True, alpha=.2, which='both')
    fig.tight_layout(); save(fig, 'v2_fig4')


# ------------------------------------------------------------ Figure 5 (new)
def figure5():
    """Plain partial product vs exact limit, B up to 1e6 (log axis)."""
    fig, ax = plt.subplots(figsize=(10, 5))
    P = sieve(10 ** 6)
    # omega via roots for speed
    def roots(p):
        # smallest primitive root
        phi = p - 1; m = phi; fac = set(); d = 2
        while d * d <= m:
            while m % d == 0: fac.add(d); m //= d
            d += 1
        if m > 1: fac.add(m)
        g = next(g for g in range(2, p) if all(pow(g, phi // q, p) != 1 for q in fac))
        z = pow(g, phi // 47, p); out = set(); y = 1
        for _ in range(46):
            y = y * z % p; out.add((1 + pow(y - 1, -1, p)) % p)
        return out
    run = {4: 0.0, 5: 0.0, 6: 0.0}; xs = []; ys = {4: [], 5: [], 6: []}
    for p in P:
        if p % 47 == 1:
            R = roots(p); om = {k: len({(r - i) % p for r in R for i in range(k)}) for k in run}
        else:
            om = {k: 0 for k in run}
        for k in run:
            run[k] += math.log1p(-om[k] / p) - k * math.log1p(-1 / p)
        if p % 47 == 1 or p < 2000 or p % 7 == 1:
            xs.append(p)
            for k in run: ys[k].append(100 * (math.exp(run[k]) / S[k] - 1))
    for k, c in [(4, '#e74c3c'), (5, '#9b59b6'), (6, '#3498db')]:
        ax.semilogx(xs, ys[k], color=c, lw=1.4, label=r'$k=%d$' % k)
    ax.axhline(0, color='k', lw=0.8); ax.axvline(1e4, color='grey', ls='--', alpha=0.7)
    ax.text(1.1e4, 45, 'v1 truncation\n$B=10^4$', fontsize=9, color='grey')
    ax.set_xlim(100, 1e6); ax.set_ylim(-60, 80)
    ax.set_xlabel('Sieve bound $B$'); ax.set_ylabel(r'$\mathfrak{S}_k(B)\,/\,\mathfrak{S}_k - 1$  (%)')
    ax.set_title('Plain partial Euler product relative to the exact (L-function) value', fontweight='bold')
    ax.legend(fontsize=10); ax.grid(True, alpha=0.2, which='both')
    save(fig, 'v2_fig5')


if __name__ == '__main__':
    print('Generating v2 figures ...')
    figure1(); figure2(); figure3(); figure4(); figure5()
    print('Done -> figures/v2/')
