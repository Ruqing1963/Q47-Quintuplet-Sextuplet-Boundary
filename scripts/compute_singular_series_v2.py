#!/usr/bin/env python3
"""
compute_singular_series_v2.py -- Exact Bateman-Horn singular series for the
k consecutive shifts of Q(n) = n^47 - (n-1)^47, k = 1..6  (Paper III, v2).

Two evaluations are produced:

(1) PLAIN partial Euler product (what v1 of the paper tabulated)
        S_k(B) = prod_{p <= B} (1 - omega_k(p)/p) / (1 - 1/p)^k .
    This converges only conditionally: it oscillates by tens of percent
    for B in [283, 10^4] and cannot by itself certify convergence.

(2) ACCELERATED (absolutely convergent) evaluation
        S_k = [ prod_{chi != chi0 (mod 47)} L(1,chi) ]^(-k) * prod_p B_k(p),
        B_k(p) = (1 - omega_k(p)/p) * (1 - p^{-f})^{-46k/f},   f = ord_47(p), p != 47,
        B_k(47) = (1 - 1/47)^{-k}.
    Because omega_k(p) = 46k * [p = 1 mod 47] for all but finitely many p,
    B_k(p) = 1 + O(p^-2) and the tail beyond B is O(k^2 / (B log B)).
    The product of L-values is computed exactly from the finite formulas for
    primitive characters modulo a prime (mpmath, 30 digits).

omega_k(p) is obtained from the 46 explicit roots  r_j = 1 + (zeta^j - 1)^{-1},
zeta a primitive 47th root of unity mod p (cross-checked by brute force).

Usage:
    python compute_singular_series_v2.py [--bound 1000000] [--write-csv]

Author: Ruqing Chen, GUT Geoservice Inc.   (v2: October 2026)
"""
import argparse, math, os, sys, time, csv
import numpy as np
from mpmath import mp, mpc, pi, exp, log, fabs, sqrt

mp.dps = 30
Q_EXP = 47
KS = [1, 2, 3, 4, 5, 6]
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')


# ----------------------------------------------------------------------------
def primes_upto(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]


def primitive_root(p):
    phi = p - 1; m = phi; fac = set(); d = 2
    while d * d <= m:
        while m % d == 0:
            fac.add(d); m //= d
        d += 1
    if m > 1:
        fac.add(m)
    for g in range(2, p):
        if all(pow(g, phi // f, p) != 1 for f in fac):
            return g
    raise RuntimeError(p)


def roots_mod_p(p):
    """The 46 roots of Q(n) mod a resonant prime p (p = 1 mod 47)."""
    zeta = pow(primitive_root(p), (p - 1) // Q_EXP, p)
    roots, z = set(), 1
    for _ in range(1, Q_EXP):
        z = z * zeta % p
        roots.add((1 + pow(z - 1, -1, p)) % p)
    assert len(roots) == Q_EXP - 1
    return roots


def omega_from_roots(roots, p, k):
    """omega_k(p) = #{r mod p : some Q(r+i), 0<=i<k, vanishes} = |U (roots - i)|."""
    return len({(r - i) % p for r in roots for i in range(k)})


def omega_bruteforce(p, k):
    Qm = lambda n: (pow(n, Q_EXP, p) - pow(n - 1, Q_EXP, p)) % p
    return p - sum(1 for r in range(p) if all(Qm(r + i) != 0 for i in range(k)))


def ord_mod_q(r):
    f, x = 1, r % Q_EXP
    while x != 1:
        x = x * r % Q_EXP; f += 1
    return f


# ----------------------------------------------------------------------------
def L1_values():
    """L(1,chi) for the 45 non-principal characters mod 47 (all primitive).
    chi_m(g^a) = exp(2 pi i m a / 46).  Even chi: m even."""
    q = Q_EXP
    g = primitive_root(q)
    ind, x = {}, 1
    for a in range(q - 1):
        ind[x] = a; x = x * g % q
    vals = []
    for m in range(1, q - 1):
        chi = lambda a, m=m: exp(2j * pi * m * ind[a % q] / (q - 1)) if a % q else mpc(0)
        tau = sum(chi(a) * exp(2j * pi * a / q) for a in range(1, q))          # Gauss sum
        if m % 2 == 0:   # even character
            L = -(tau / q) * sum(chi(a).conjugate() * log(fabs(1 - exp(2j * pi * a / q)))
                                 for a in range(1, q))
        else:            # odd character
            L = (pi * 1j * tau / q ** 2) * sum(chi(a).conjugate() * a for a in range(1, q))
        vals.append(L)
    # internal check: the quadratic character (m = 23) gives L(1,chi) = pi h(-47)/sqrt(47), h = 5
    quad = vals[22]
    assert abs(quad - pi * 5 / sqrt(47)) < mp.mpf(10) ** -20, quad
    prod = mpc(1)
    for v in vals:
        prod *= v
    assert abs(prod.imag) < 1e-20
    return float(prod.real), vals


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bound', '-B', type=int, default=10 ** 6)
    ap.add_argument('--write-csv', action='store_true', help='write data/*.csv (v2 files)')
    args = ap.parse_args()
    B = args.bound
    t0 = time.time()

    print("Cross-check omega_k(p) (roots vs brute force) at the first resonant primes:")
    for p in (283, 659, 941, 1129, 1223):
        R = roots_mod_p(p)
        a = [omega_from_roots(R, p, k) for k in KS]
        b = [omega_bruteforce(p, k) for k in KS]
        assert a == b, (p, a, b)
        print(f"  p={p:>5}: omega_1..6 = {a}")

    Lprod, _ = L1_values()
    print(f"\nprod_(chi != chi0) L(1,chi) = {Lprod:.12f}")

    primes = primes_upto(B)
    print(f"primes up to {B}: {len(primes)}")

    logS = {k: 0.0 for k in KS}        # plain partial product
    logA = {k: 0.0 for k in KS}        # accelerated product of B_k(p)
    plain_checkpoints = [100, 282, 283, 1000, 5000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
    accel_checkpoints = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
    plain_rows, accel_rows, omega_rows = [], [], []
    n_res = 0
    pc = ac = 0

    def flush_checkpoints(p_done):
        nonlocal pc, ac
        while pc < len(plain_checkpoints) and plain_checkpoints[pc] <= p_done:
            plain_rows.append((plain_checkpoints[pc], {k: math.exp(logS[k]) for k in KS}))
            pc += 1
        while ac < len(accel_checkpoints) and accel_checkpoints[ac] <= p_done:
            accel_rows.append((accel_checkpoints[ac], {k: math.exp(logA[k]) * Lprod ** (-k) for k in KS}))
            ac += 1

    prev = 1
    for p in primes:
        # checkpoints lying strictly between the previous prime and this one
        flush_checkpoints(p - 1)
        if p == Q_EXP:
            for k in KS:
                logS[k] += -k * math.log1p(-1 / p)
                logA[k] += -k * math.log1p(-1 / p)
        else:
            f = ord_mod_q(p)
            if f == 1:
                n_res += 1
                om = {k: omega_from_roots(roots_mod_p(p), p, k) for k in KS}
                if p < 10 ** 4:
                    omega_rows.append((p, om))
            else:
                om = {k: 0 for k in KS}
            for k in KS:
                a = math.log1p(-om[k] / p)
                logS[k] += a - k * math.log1p(-1 / p)
                logA[k] += a - (46 * k / f) * math.log1p(-float(p) ** (-f))
        prev = p
    flush_checkpoints(B)

    S = {k: math.exp(logA[k]) * Lprod ** (-k) for k in KS}

    print(f"\nresonant primes (p = 1 mod 47) up to {B}: {n_res}  "
          f"(Dirichlet expectation {len(primes)/46:.0f})")

    print("\nomega_k(p) at resonant primes below 10^4 (k = 1..6):")
    for p, om in omega_rows:
        flag = '' if om[6] == 276 else '   <- overlapping shifted roots'
        print(f"  p={p:>5}: {[om[k] for k in KS]}{flag}")

    print("\nPLAIN partial Euler product S_k(B):")
    print(f"{'B':>9} " + " ".join(f"{'S'+str(k):>12}" for k in KS))
    for b, row in plain_rows:
        print(f"{b:>9} " + " ".join(f"{row[k]:>12.6g}" for k in KS))

    print("\nACCELERATED evaluation (converges absolutely):")
    print(f"{'B':>9} " + " ".join(f"{'S'+str(k):>12}" for k in KS))
    for b, row in accel_rows:
        print(f"{b:>9} " + " ".join(f"{row[k]:>12.6g}" for k in KS))

    print("\n=== EXACT SINGULAR SERIES (accelerated, B = %d) ===" % B)
    for k in KS:
        print(f"  S_{k} = {S[k]:.6g}")
    print("  ratios S_{k+1}/S_k:",
          "  ".join(f"{k}->{k+1}: {S[k+1]/S[k]:.4f}" for k in KS[:-1]))
    print(f"  S_6/S_4 = {S[6]/S[4]:.3f}")
    print(f"elapsed {time.time()-t0:.1f}s")

    if args.write_csv:
        os.makedirs(DATA, exist_ok=True)
        with open(os.path.join(DATA, 'singular_series_exact.csv'), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['k', 'S_k', 'S_k_over_S_km1', 'method', 'bound'])
            for k in KS:
                w.writerow([k, f"{S[k]:.6g}", f"{S[k]/S[k-1]:.5f}" if k > 1 else '',
                            'L-function accelerated Euler product', B])
        with open(os.path.join(DATA, 'singular_series_convergence_v2.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sieve_bound', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'note'])
            notes = {100: 'no resonant prime yet', 282: 'just before first cliff',
                     283: 'first cliff (p=283)', 1000: 'after partial recovery',
                     5000: 'oscillating', 10 ** 4: 'v1 truncation point',
                     10 ** 5: '', 10 ** 6: '', 10 ** 7: ''}
            for b, row in plain_rows:
                w.writerow([b] + [f"{row[k]:.6g}" for k in KS] + [notes.get(b, '')])
            w.writerow(['exact'] + [f"{S[k]:.6g}" for k in KS] + ['accelerated limit'])
        with open(os.path.join(DATA, 'accelerated_convergence.csv'), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['bound', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
            for b, row in accel_rows:
                w.writerow([b] + [f"{row[k]:.7g}" for k in KS])
        with open(os.path.join(DATA, 'resonant_primes_omega_v2.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['p', 'omega_1', 'omega_2', 'omega_3', 'omega_4', 'omega_5', 'omega_6',
                        'sigma_4', 'sigma_5', 'sigma_6'])
            for p, om in omega_rows:
                sig = [(1 - om[k] / p) / (1 - 1 / p) ** k for k in (4, 5, 6)]
                w.writerow([p] + [om[k] for k in KS] + [f"{s:.6f}" for s in sig])
        with open(os.path.join(DATA, 'L1_product.txt'), 'w') as f:
            f.write(f"prod_(chi != chi0 mod 47) L(1,chi) = {Lprod:.15f}\n")
        print("CSV files written to data/ (singular_series_exact.csv, singular_series_convergence_v2.csv, "
              "accelerated_convergence.csv, resonant_primes_omega_v2.csv, L1_product.txt)")


if __name__ == '__main__':
    main()
