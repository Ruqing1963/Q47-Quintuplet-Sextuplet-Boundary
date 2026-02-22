#!/usr/bin/env python3
"""
compute_singular_series.py — Compute Bateman-Horn singular series for
Q(n) = n^47 - (n-1)^47 with CORRECT root counts.

Key finding: Q(n) has a bifurcated root structure:
  - ω₁(p) = 0  for all primes p ≢ 1 (mod 47)
  - ω₁(p) = 46 for all primes p ≡ 1 (mod 47)  ("resonant primes")

This script computes 𝔖_k for k = 4, 5, 6 via Euler products to
a given sieve bound B, incorporating the correct ω_k(p) at each
resonant prime.

Usage:
    python compute_singular_series.py [--bound 10000]

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import argparse
import sys


def sieve_primes(N):
    """Return list of primes up to N via Sieve of Eratosthenes."""
    s = [True] * (N + 1)
    s[0] = s[1] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]:
            for j in range(i * i, N + 1, i):
                s[j] = False
    return [i for i in range(2, N + 1) if s[i]]


def Q_mod(n, p):
    """Compute Q(n) mod p = (n^47 - (n-1)^47) mod p."""
    return (pow(n, 47, p) - pow(n - 1, 47, p)) % p


def compute_omega_k(p, k):
    """
    Compute ω_k(p): number of residues r in [0, p) where at least one
    of Q(r), Q(r+1), ..., Q(r+k-1) vanishes mod p.
    """
    admissible = 0
    for r in range(p):
        if all(Q_mod(r + i, p) != 0 for i in range(k)):
            admissible += 1
    return p - admissible


def main():
    parser = argparse.ArgumentParser(
        description='Compute Bateman-Horn singular series for Q(n)=n^47-(n-1)^47')
    parser.add_argument('--bound', '-B', type=int, default=10000,
                        help='Sieve bound (default: 10000)')
    args = parser.parse_args()

    B = args.bound
    primes = sieve_primes(B)
    print(f"Computing singular series to sieve bound B = {B}")
    print(f"  {len(primes)} primes to process\n")

    # Compute Euler products
    S = {4: 1.0, 5: 1.0, 6: 1.0}
    resonant_count = 0

    print(f"{'p':>6} {'p%47':>4} {'ω₁':>4} {'ω₄':>4} {'ω₅':>4} {'ω₆':>4} "
          f"{'σ₄':>8} {'σ₅':>8} {'σ₆':>8}  {'S₄':>10} {'S₅':>10} {'S₆':>12}")
    print("-" * 100)

    milestones = {100, 200, 283, 500, 1000, 2000, 5000, B}

    for p in primes:
        # Single-polynomial roots
        omega1 = sum(1 for r in range(p) if Q_mod(r, p) == 0)

        # k-constellation roots
        omegas = {}
        for k in [4, 5, 6]:
            omegas[k] = compute_omega_k(p, k)

        # Local factors
        for k in [4, 5, 6]:
            factor = (1 - omegas[k] / p) / (1 - 1 / p) ** k
            S[k] *= factor

        is_resonant = (p % 47 == 1)
        if is_resonant:
            resonant_count += 1

        # Print resonant primes and milestones
        if is_resonant or p in milestones:
            sigma4 = (1 - omegas[4] / p) / (1 - 1 / p) ** 4
            sigma5 = (1 - omegas[5] / p) / (1 - 1 / p) ** 5
            sigma6 = (1 - omegas[6] / p) / (1 - 1 / p) ** 6
            tag = " *** RESONANT" if is_resonant else ""
            print(f"{p:>6} {p % 47:>4} {omega1:>4} {omegas[4]:>4} "
                  f"{omegas[5]:>4} {omegas[6]:>4} "
                  f"{sigma4:>8.4f} {sigma5:>8.4f} {sigma6:>8.4f}  "
                  f"{S[4]:>10.1f} {S[5]:>10.1f} {S[6]:>12.1f}{tag}")

    print("-" * 100)
    print(f"\n=== RESULTS (B = {B}) ===")
    print(f"  Resonant primes (p ≡ 1 mod 47): {resonant_count}")
    print(f"  Non-resonant primes: {len(primes) - resonant_count}")
    print(f"\n  𝔖₄(B) = {S[4]:.2f}")
    print(f"  𝔖₅(B) = {S[5]:.2f}")
    print(f"  𝔖₆(B) = {S[6]:.2f}")
    print(f"\n  𝔖₅/𝔖₄ = {S[5] / S[4]:.4f}")
    print(f"  𝔖₆/𝔖₅ = {S[6] / S[5]:.4f}")
    print(f"  𝔖₆/𝔖₄ = {S[6] / S[4]:.4f}")

    # Calibrate against empirical S4 = 6385
    cal = 6385 / S[4]
    print(f"\n  Calibration factor (6385/S₄): {cal:.6f}")
    print(f"  𝔖₄ (calibrated) = {6385}")
    print(f"  𝔖₅ (calibrated) = {S[5] * cal:.0f}")
    print(f"  𝔖₆ (calibrated) = {S[6] * cal:.0f}")


if __name__ == '__main__':
    main()
