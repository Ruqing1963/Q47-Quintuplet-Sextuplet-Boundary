#!/usr/bin/env python3
"""
generate_figures.py — Reproduce all 4 figures for Paper III
(Prime Quintuplets and the Sextuplet Boundary).

Requires: matplotlib, numpy, scipy
Data:     ../data/ CSV files

Usage:
    python generate_figures.py

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from collections import Counter

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11, 'axes.labelsize': 12,
    'axes.titlesize': 13, 'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'mathtext.fontset': 'cm',
})

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')
OUT  = os.path.join(os.path.dirname(__file__), '..', 'paper', 'figures')
os.makedirs(OUT, exist_ok=True)

# Corrected singular series
S4, S5, S6 = 6385, 57108, 519756

# Load data
quints = np.array([35676017721, 64482563907, 73292417435, 116255850744,
                    147743683226, 159430471996, 182501065420])

# Try to load quadruplet data from Part II repo or use inline
quad_file = os.path.join(DATA, '..', '..', 'Q47-Deep-Space-Quadruplet-Census',
                         'data', 'quadruplets_742.csv')
if os.path.exists(quad_file):
    quad_starts = []
    with open(quad_file) as f:
        next(f)
        for line in f:
            quad_starts.append(int(line.strip().split(',')[1]))
    quad_starts = np.array(sorted(quad_starts))
else:
    # Fallback: generate synthetic positions matching density profile
    print("  Note: quadruplet data not found, using quintuplet-only mode")
    quad_starts = None


def Q_mod(n, p):
    return (pow(n, 47, p) - pow(n - 1, 47, p)) % p


def sieve_primes(N):
    s = [True] * (N + 1)
    s[0] = s[1] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]:
            for j in range(i * i, N + 1, i):
                s[j] = False
    return [i for i in range(2, N + 1) if s[i]]


# ── Figure 1: Quintuplet field ─────────────────────────────────
def figure1():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6),
        height_ratios=[2, 1.5], gridspec_kw={'hspace': 0.25})

    if quad_starts is not None:
        ax1.scatter(quad_starts / 1e9, np.ones_like(quad_starts),
                    c='#AAAAAA', s=6, alpha=0.4, marker='|', linewidths=0.8)

    for i, q in enumerate(quints):
        ax1.scatter(q / 1e9, 1.0, c='gold', s=250, marker='*',
                    edgecolors='darkorange', linewidths=1.0, zorder=5)
        ax1.annotate(f'Q{i+1}', xy=(q / 1e9, 1.0), xytext=(q / 1e9, 1.25),
                     fontsize=8, ha='center', color='#8B4513', fontweight='bold')

    ax1.set_xlim(-5, 210); ax1.set_ylim(0.5, 1.6); ax1.set_yticks([])
    ax1.set_xlabel(r'$n$ (billions)')
    ax1.set_title('Figure 1.  Quintuplet positions within the quadruplet field',
                  fontweight='bold', pad=10)

    # Cumulative count
    N_q = np.sort(quints)
    C_q = np.arange(1, len(N_q) + 1)
    ax2.step(np.concatenate([[0], N_q / 1e9, [200]]),
             np.concatenate([[0], C_q, [7]]),
             color='darkorange', linewidth=2.5, where='post', label='Observed')

    N_fit = np.linspace(1e9, 2.1e11, 500)
    C_bh = [S5 * integrate.quad(lambda t: 1.0 / (46 * np.log(t))**5, 2, n)[0]
            for n in N_fit]
    ax2.plot(N_fit / 1e9, C_bh, '--', color='#27ae60', linewidth=1.5,
             label=r'Bateman-Horn ($\mathfrak{S}_5 \approx 57{,}100$)')

    ax2.set_xlabel(r'$N$ (billions)'); ax2.set_ylabel(r'$C_5(N)$')
    ax2.set_xlim(0, 210); ax2.set_ylim(0, 9)
    ax2.legend(fontsize=10, loc='upper left'); ax2.grid(True, alpha=0.2)
    plt.savefig(os.path.join(OUT, 'p2_fig1_v2.png'), dpi=300); plt.close()
    print("  Figure 1 saved.")


# ── Figure 2: Staircase ────────────────────────────────────────
def figure2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    labels = ['Solitary\n(k=1)', 'Pair\n(k=2)', 'Triplet\n(k=3)',
              'Quad\n(k=4)', 'Quint\n(k=5)', 'Sext\n(k=6)']
    counts = [18121562, 173351, 1749, 742, 7, 0]
    log_c = [np.log10(max(c, 0.047)) for c in counts]
    colors = ['#3498db', '#2ecc71', '#e67e22', '#e74c3c', '#9b59b6', '#95a5a6']
    bars = ax1.bar(range(6), log_c, color=colors, alpha=0.8,
                   edgecolor='#333', linewidth=0.8)
    for i, (bar, c) in enumerate(zip(bars, counts)):
        lbl = f'{c:,}' if c > 0 else f'0\n(E=0.047)'
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                 lbl, ha='center', fontsize=8, fontweight='bold')
    ax1.set_xticks(range(6)); ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel(r'$\log_{10}$ (count)')
    ax1.set_title('(a) Morphological hierarchy', fontweight='bold')
    ax1.set_ylim(-2, 8); ax1.grid(True, alpha=0.15, axis='y')

    n_vals = np.linspace(1e10, 2e11, 300)
    p5 = (S5 / S4) / (46 * np.log(n_vals))
    p6 = (S6 / S5) / (46 * np.log(n_vals))
    ax2.plot(n_vals / 1e9, p5 * 100, '-', color='#9b59b6', linewidth=2,
             label=r'$P(\mathrm{quint}|\mathrm{quad})$')
    ax2.plot(n_vals / 1e9, p6 * 100, '-', color='#e74c3c', linewidth=2,
             label=r'$P(\mathrm{sext}|\mathrm{quint})$')
    for q in quints:
        pq = (S5 / S4) / (46 * np.log(q))
        ax2.plot(q / 1e9, pq * 100, '*', color='gold', markersize=12,
                 markeredgecolor='darkorange', markeredgewidth=0.8, zorder=5)
    ax2.set_xlabel(r'$n$ (billions)')
    ax2.set_ylabel('Conditional probability (%)')
    ax2.set_title('(b) Extension probabilities', fontweight='bold')
    ax2.legend(fontsize=10); ax2.grid(True, alpha=0.2); ax2.set_xlim(0, 205)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p2_fig2_v2.png'), dpi=300); plt.close()
    print("  Figure 2 saved.")


# ── Figure 3: Sextuplet prediction ─────────────────────────────
def figure3():
    fig, ax = plt.subplots(figsize=(10, 6))
    N_range = np.logspace(10, 14, 500)
    E5v = [S5 * integrate.quad(lambda t: 1.0 / (46 * np.log(t))**5, 2, N)[0]
           for N in N_range]
    E6v = [S6 * integrate.quad(lambda t: 1.0 / (46 * np.log(t))**6, 2, N)[0]
           for N in N_range]

    ax.loglog(N_range, E6v, '-', color='#e74c3c', linewidth=2.5,
              label=r'$E[C_6(N)]$ ($\mathfrak{S}_6 \approx 520{,}000$)')
    ax.loglog(N_range, E5v, '-', color='#9b59b6', linewidth=2,
              label=r'$E[C_5(N)]$ ($\mathfrak{S}_5 \approx 57{,}100$)')

    ax.axvline(x=2e11, color='#3498db', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.text(2.3e11, 0.005, 'Current limit\n$N=2\\times 10^{11}$',
            fontsize=9, color='#3498db')
    ax.axhline(y=1, color='grey', linestyle=':', alpha=0.6)
    ax.text(1.2e10, 1.3, r'$E[C_6]=1$', fontsize=10, color='grey')

    ax.plot(2e11, 7, 'o', color='#9b59b6', markersize=10,
            markeredgecolor='black', markeredgewidth=1, zorder=5,
            label='Observed: 7 quintuplets')
    ax.plot(2e11, 0.047, 's', color='#e74c3c', markersize=10,
            markeredgecolor='black', markeredgewidth=1, zorder=5,
            label='Observed: 0 sextuplets')

    # Find crossing
    def f(logN):
        N = 10**logN
        I6 = integrate.quad(lambda t: 1.0 / (46 * np.log(t))**6, 2, N)[0]
        return S6 * I6 - 1
    logN_cross = optimize.brentq(f, 12, 15)
    N_cross = 10**logN_cross
    ax.plot(N_cross, 1, 'D', color='#e74c3c', markersize=10,
            markeredgecolor='black', markeredgewidth=1, zorder=5)
    ax.annotate(f'$N^* \\approx {N_cross:.2e}$', xy=(N_cross, 1),
                xytext=(N_cross * 2, 0.3), fontsize=10, color='#c0392b',
                arrowprops=dict(arrowstyle='->', color='#c0392b'),
                fontweight='bold')

    ax.set_xlabel(r'Search bound $N$'); ax.set_ylabel('Expected count')
    ax.set_title('Figure 3.  Predicted constellation counts (corrected)',
                 fontweight='bold', pad=10)
    ax.legend(fontsize=9.5, loc='upper left')
    ax.set_xlim(5e9, 1.5e14); ax.set_ylim(0.003, 500)
    ax.grid(True, alpha=0.2, which='both')
    plt.savefig(os.path.join(OUT, 'p2_fig3_v2.png'), dpi=300); plt.close()
    print(f"  Figure 3 saved. Sextuplet boundary: N* = {N_cross:.3e}")


# ── Figure 4: Euler product structure ──────────────────────────
def figure4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    primes_show = sieve_primes(50)
    f4, f5, f6 = [], [], []
    for p in primes_show:
        for k, fl in [(4, f4), (5, f5), (6, f6)]:
            adm = sum(1 for r in range(p)
                      if all(Q_mod(r + i, p) != 0 for i in range(k)))
            omega = p - adm
            fl.append((1 - omega / p) / (1 - 1 / p)**k)

    x = np.arange(len(primes_show))
    w = 0.25
    ax1.bar(x - w, f4, w, color='#e74c3c', alpha=0.8, label=r'$k=4$')
    ax1.bar(x, f5, w, color='#9b59b6', alpha=0.8, label=r'$k=5$')
    ax1.bar(x + w, f6, w, color='#3498db', alpha=0.8, label=r'$k=6$')
    ax1.set_xticks(x); ax1.set_xticklabels(primes_show, fontsize=8)
    ax1.set_xlabel(r'Prime $p$')
    ax1.set_ylabel(r'BH local factor $\sigma_k(p)$')
    ax1.set_title(r'(a) Local factors ($p<50$: $\omega_k=0$)', fontweight='bold')
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.2, axis='y')

    # Cumulative product with cliffs
    all_p = sieve_primes(2000)
    S4c, S5c, S6c = [], [], []
    p4, p5, p6 = 1.0, 1.0, 1.0
    for p in all_p:
        for k in [4, 5, 6]:
            adm = sum(1 for r in range(p)
                      if all(Q_mod(r + i, p) != 0 for i in range(k)))
            omega = p - adm
            fac = (1 - omega / p) / (1 - 1 / p)**k
            if k == 4: p4 *= fac
            elif k == 5: p5 *= fac
            else: p6 *= fac
        S4c.append(p4); S5c.append(p5); S6c.append(p6)

    ax2.semilogy(all_p, S4c, '-', color='#e74c3c', linewidth=2, label=r'$\mathfrak{S}_4$')
    ax2.semilogy(all_p, S5c, '-', color='#9b59b6', linewidth=2, label=r'$\mathfrak{S}_5$')
    ax2.semilogy(all_p, S6c, '-', color='#3498db', linewidth=2, label=r'$\mathfrak{S}_6$')
    ax2.axvline(x=283, color='#333', linestyle=':', alpha=0.5)
    ax2.text(300, 2e5, '$p=283$\n(first cliff)', fontsize=8, color='#333')
    ax2.axvline(x=659, color='#333', linestyle=':', alpha=0.3)
    ax2.text(680, 5e4, '659', fontsize=7, color='#666')
    ax2.set_xlabel(r'Sieve bound $B$')
    ax2.set_ylabel(r'Cumulative $\mathfrak{S}_k(B)$')
    ax2.set_title('(b) Euler product with periodic obstruction', fontweight='bold')
    ax2.legend(fontsize=10, loc='upper left'); ax2.grid(True, alpha=0.2, which='both')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p2_fig4_v2.png'), dpi=300); plt.close()
    print("  Figure 4 saved.")


if __name__ == '__main__':
    print("Generating figures for Paper III (Quintuplets & Sextuplet Boundary)...")
    figure1()
    figure2()
    figure3()
    figure4()
    print("Done.")
