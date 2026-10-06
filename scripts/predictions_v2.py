#!/usr/bin/env python3
"""
predictions_v2.py -- Bateman-Horn predictions for Paper III (v2) from the exact
singular series in data/singular_series_exact.csv.

Produces (and prints):
  * E[C_k(N)] for the pioneer zone (N = 2e9, k = 1..4) and the full range
    (N = 2e11, k = 4..6), with and without the ln 47 refinement of ln Q(t);
  * the empirical calibration S_4(emp) = 742 / I_4 and its ratio to the exact S_4;
  * inter-level suppression factors R_{k->k+1}(N) = (S_k/S_{k+1}) (I_k/I_{k+1});
  * the sextuplet boundary N* (E[C_6(N*)] = 1) and a prediction table;
  * the conditional extension model summed over the 742 quadruplets;
  * KS test of the 7 quintuplet positions against the BH cumulative;
  * corrected digit counts of the quintuplets.

Writes: data/hierarchy_validation.csv, data/sextuplet_predictions_v2.csv,
        data/suppression_ratios_v2.csv, data/quintuplets_7_v2.csv

Usage:  python predictions_v2.py
"""
import csv, math, os
import numpy as np
from scipy import integrate, optimize
from scipy.stats import kstest, poisson

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
LN47 = math.log(47)

# ---------------------------------------------------------------- inputs
S = {}
with open(os.path.join(DATA, 'singular_series_exact.csv')) as f:
    for row in csv.DictReader(f):
        S[int(row['k'])] = float(row['S_k'])

quads = []
with open(os.path.join(DATA, 'quadruplets_742.csv')) as f:
    for row in csv.DictReader(f):
        quads.append(int(row['starting_n']))
quads = np.array(sorted(quads))
assert len(quads) == 742

quints = np.array([35676017721, 64482563907, 73292417435, 116255850744,
                   147743683226, 159430471996, 182501065420])

# Observed counts of MAXIMAL runs (Parts I and II): a quintuplet is listed once,
# not also as a quadruplet.  Bateman-Horn counts POSITIONS n with Q(n..n+k-1)
# all prime, so a maximal run of length m contributes (m-k+1) positions.
# Part I categories are read as runs of EXACT length m (solitary = isolated prime,
# pair = exactly two, ...).  Part II's 742 quadruplets include the 7 quintuplets,
# i.e. 735 runs of exact length 4 and 7 of exact length 5.
RUNS = {2e9: {1: 18121562, 2: 173351, 3: 1749, 4: 15, 5: 0},
        2e11: {4: 735, 5: 7, 6: 0}}
def positions(N, k):
    return sum((m - k + 1) * c for m, c in RUNS[N].items() if m >= k)
def catalogue(N, k):      # runs of length >= k (the number a census reports)
    return sum(c for m, c in RUNS[N].items() if m >= k)
OBS = {(N, k): positions(N, k) for N in RUNS for k in RUNS[N]}
OBS_RUNS = {(N, k): catalogue(N, k) for N in RUNS for k in RUNS[N]}


# ---------------------------------------------------------------- BH integrals
def lnQ(t):
    return 46 * math.log(t) + LN47          # ln Q(t) up to O(1/t)


def I(k, N, refined=True):
    g = (lambda t: 1.0 / lnQ(t) ** k) if refined else (lambda t: 1.0 / (46 * math.log(t)) ** k)
    return integrate.quad(g, 2, N, limit=500)[0]


def E(k, N, refined=True, Sk=None):
    return (S[k] if Sk is None else Sk) * I(k, N, refined)


print("Exact singular series:", {k: round(v, 1) for k, v in S.items()})

# ---------------------------------------------------------------- hierarchy validation
print("\n=== Hierarchy validation: Bateman-Horn (no free parameters) vs observed ===")
rows = []
for (N, k), obs in OBS.items():
    e_ref, e_plain = E(k, N), E(k, N, refined=False)
    sd = math.sqrt(e_ref)
    z = (obs - e_ref) / sd
    runs = OBS_RUNS[(N, k)]
    rows.append([f"{N:.0e}", k, runs, obs, f"{e_ref:.4g}", f"{e_plain:.4g}", f"{z:+.2f}",
                 f"{100 * (obs / e_ref - 1):+.2f}"])
    print(f"  N={N:.0e} k={k}: catalogue (runs>=k) {runs:>10,}  positions {obs:>10,}   "
          f"E(refined) {e_ref:12,.3f}   E(46 ln t) {e_plain:12,.3f}   "
          f"dev {100 * (obs / e_ref - 1):+.2f}%  z = {z:+.2f}")
with open(os.path.join(DATA, 'hierarchy_validation.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'k', 'observed_catalogue_runs_ge_k', 'observed_positions', 'E_BH_refined',
                'E_BH_46lnt', 'z_score', 'deviation_percent'])
    w.writerows(rows)

I4 = I(4, 2e11)
P4 = OBS[(2e11, 4)]
S4_emp = P4 / I4
print(f"\nEmpirical calibration S_4(emp) = {P4} / I_4 = {S4_emp:.1f}  "
      f"(v1 convention 742/(46 ln t): {742 / I(4, 2e11, False):.1f});  exact S_4 = {S[4]:.1f};  "
      f"ratio = {S4_emp / S[4]:.4f}  (Poisson 1 sigma = {100 / math.sqrt(P4):.1f} %)")

# ---------------------------------------------------------------- suppression factors
print("\n=== Suppression factors R_{k->k+1}(N) = (S_k/S_{k+1}) * (I_k/I_{k+1}) ===")
emp = {(2e9, 1): OBS[(2e9, 1)] / OBS[(2e9, 2)], (2e9, 2): OBS[(2e9, 2)] / OBS[(2e9, 3)],
       (2e9, 3): OBS[(2e9, 3)] / OBS[(2e9, 4)], (2e11, 4): OBS[(2e11, 4)] / OBS[(2e11, 5)]}
srows = []
for N in (2e9, 2e11):
    for k in range(1, 6):
        R = (S[k] / S[k + 1]) * (I(k, N) / I(k + 1, N))
        e = emp.get((N, k))
        srows.append([f"{N:.0e}", f"{k}->{k+1}", f"{R:.1f}", f"{e:.1f}" if e else ''])
        print(f"  N={N:.0e}  {k}->{k+1}: predicted {R:6.1f}   observed {e if e is None else round(e,1)}")
with open(os.path.join(DATA, 'suppression_ratios_v2.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'transition', 'R_predicted', 'R_observed']); w.writerows(srows)

# ---------------------------------------------------------------- sextuplet boundary
print("\n=== Sextuplet boundary ===")
def boundary(S6):
    return 10 ** optimize.brentq(lambda lg: S6 * I(6, 10 ** lg) - 1, 12, 15)
Nstar = boundary(S[6])
Nstar_cal = boundary(S[6] * S4_emp / S[4])
Nstar_plain = 10 ** optimize.brentq(lambda lg: S[6] * I(6, 10 ** lg, False) - 1, 12, 15)
digits = len(str(int(Nstar) ** 47 - (int(Nstar) - 1) ** 47))
print(f"  exact S_6:                N* = {Nstar:.3e}  ({Nstar/2e11:.1f} x current range), Q has {digits} digits")
print(f"  calibrated S_6 (x{S4_emp/S[4]:.3f}): N* = {Nstar_cal:.3e}")
print(f"  46 ln t convention:       N* = {Nstar_plain:.3e}")
print(f"  P(>=1 sextuplet by 2e11) = {1 - math.exp(-E(6, 2e11)):.3%}")
prows = []
for N in (2e11, 1e12, 2e12, 5e12, 1e13, Nstar, 2e13, 5e13, 1e14):
    e5, e6 = E(5, N), E(6, N)
    prows.append([f"{N:.3e}", f"{e5:.1f}", f"{e6:.3f}", f"{1 - math.exp(-e6):.3f}"])
    print(f"  N={N:.3e}: E[C5]={e5:8.1f}  E[C6]={e6:.3f}  P(>=1 sext)={1 - math.exp(-e6):.3f}")
with open(os.path.join(DATA, 'sextuplet_predictions_v2.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'E_C5', 'E_C6', 'P_at_least_one_sextuplet']); w.writerows(prows)

# ---------------------------------------------------------------- conditional model
print("\n=== Conditional extension model over the 742 quadruplets ===")
ratio = S[5] / S[4]
mean_ln = float(np.mean(np.log(quads)))
P_each = ratio / np.array([lnQ(n + 4) for n in quads])
E_cond = float(P_each.sum())
E_cond_mean = 742 * ratio / (46 * mean_ln)
print(f"  S5/S4 = {ratio:.4f};  <ln n> over quadruplets = {mean_ln:.3f}")
print(f"  E_cond (sum of per-quadruplet probabilities, ln Q refined) = {E_cond:.3f}")
print(f"  E_cond (742 * ratio / (46 <ln n>))                       = {E_cond_mean:.3f}")
print(f"  naive independence (ratio = 1): {E_cond / ratio:.3f}")
E5 = E(5, 2e11)
print(f"  direct E[C5] = {E5:.3f};  P(X >= 7 | {E5:.2f}) = {1 - poisson.cdf(6, E5):.3f};  "
      f"1-sigma band [{E5 - math.sqrt(E5):.2f}, {E5 + math.sqrt(E5):.2f}]")

# ---------------------------------------------------------------- KS test
print("\n=== KS test of quintuplet positions ===")
I5N = I(5, 2e11)
cdf = np.vectorize(lambda x: I(5, x) / I5N if x > 2 else 0.0)
ks_bh = kstest(quints.astype(float), cdf)
ks_u = kstest(quints / 2e11, 'uniform')
print(f"  vs Bateman-Horn cumulative: D = {ks_bh.statistic:.3f}, p = {ks_bh.pvalue:.2f}")
print(f"  vs uniform on [0, 2e11]:    D = {ks_u.statistic:.3f}, p = {ks_u.pvalue:.2f}")
gaps = np.diff(quints) / 1e9
print(f"  gaps (1e9): {np.round(gaps, 1)}  mean {gaps.mean():.1f}  sd(pop) {gaps.std():.1f}  sd(sample) {gaps.std(ddof=1):.1f}")

# ---------------------------------------------------------------- quintuplet table
print("\n=== Quintuplet catalogue with exact digit counts ===")
ranks = {n: i + 1 for i, n in enumerate(quads)}
with open(os.path.join(DATA, 'quintuplets_7_v2.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['index', 'starting_n', 'digits_Q_n', 'digits_Q_n_plus_4', 'n_over_1e11',
                'quad_rank', 'gap_from_previous_1e9'])
    for i, n in enumerate(quints):
        d0 = len(str(int(n) ** 47 - (int(n) - 1) ** 47))
        d4 = len(str((int(n) + 4) ** 47 - (int(n) + 3) ** 47))
        gap = '' if i == 0 else f"{(n - quints[i-1]) / 1e9:.1f}"
        w.writerow([i + 1, int(n), d0, d4, f"{n / 1e11:.4f}", ranks.get(int(n), ''), gap])
        print(f"  #{i+1} n={int(n):>15,}  digits {d0}/{d4}  quad rank #{ranks.get(int(n))}  gap {gap}")
print("\nCSV files written.")
