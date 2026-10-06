# Changelog — Version 2 (October 2026)

**Version 2 DOI: [10.5281/zenodo.23178170](https://doi.org/10.5281/zenodo.23178170)**
(Version 1 DOI: [10.5281/zenodo.18728917](https://doi.org/10.5281/zenodo.18728917))

Version 1 of *Prime Quintuplets and the Sextuplet Boundary for Q(n) = n⁴⁷ − (n−1)⁴⁷*
(Zenodo 10.5281/zenodo.18728917, February 2026) was re-audited line by line against
its own scripts and against independent computations. The headline conclusions
survive (7 quintuplets and 0 sextuplets are what Bateman–Horn predicts; the
sextuplet boundary is near 10¹³), but several numbers and one qualitative claim
were wrong. This file lists every change.

## Errors corrected

| # | Location (v1) | v1 statement | Correct | Severity |
|---|---|---|---|---|
| 1 | Table 2 (partial singular series) | 𝔖₄(B) = 3,389 / 3,601 / **2,026** / 3,817 / 5,694 for B = 100 / 283⁻ / 283 / 1000 / 5000 (and the 𝔖₅, 𝔖₆ columns) | 4,772 / 10,780 / **6,066** / 7,947 / 6,668 (𝔖₅, 𝔖₆ likewise). Only the B = 10,000 row (6,392 / 57,173 / 520,348) was the output of `compute_singular_series.py`. Fig. 4(b) of v1 already showed the correct curve and contradicted the table. | **Major** |
| 2 | Table 2, "root-free" column; §2.3; Eq. (3) discussion | 𝔖₄^free(1000) = 8,604; "root-free" ratios 9.6, 9.8, 94.4; "35 % reduction" | The root-free product ∏(1−1/p)⁻ᵏ diverges like (e^γ ln B)ᵏ (23,270 at B = 1000; 72,773 at B = 10⁴); its ratios 𝔖₅/𝔖₄ and 𝔖₆/𝔖₅ are identical by construction. The comparison is withdrawn. | Major |
| 3 | §3.1 | "Calibration factor 6385/6392 = 0.999, indicating near-complete convergence by B = 10,000" | The plain product oscillates by ±20 % up to B ≈ 10³ and ±1 % at B ≈ 10⁵; its value at 10⁴ is 2–3 % *below* the limit. The agreement was coincidental. v2 evaluates 𝔖ₖ exactly via Dirichlet L-functions: 𝔖₄ = 6,514, 𝔖₅ = 58,550, 𝔖₆ = 535,500 (±10⁻⁴). | **Major** |
| 4 | Abstract, §1, §7, Conclusion | "stable suppression factor ≈ 127 per level … independent of the polynomial's root distribution"; Table 6 mean "≈ 110" | R_{k→k+1}(N) = (I_k/I_{k+1}) / (𝔖_{k+1}/𝔖_k) ≈ (46 ln N + ln 47)/9 grows logarithmically: ≈ 103 at N = 2×10⁹ (observed 104.4, 99.4, 118.6) and ≈ 127 at N = 2×10¹¹ (observed 107 from 749/7). | Major (qualitative) |
| 5 | Table 3, Abstract, §6 | Quintuplet digits 487–519; first sextuplet "∼605 digits" | 488–520 (v1 and the Part II catalogue used ⌊log₁₀ Q⌋ instead of ⌊log₁₀ Q⌋+1); Q(N*) has 601 digits. | Minor |
| 6 | Table 3 | Column "Gap to next" | Values are gaps **from the previous** quintuplet (the CSV `quintuplets_7.csv` was shifted by one row relative to the table). | Minor |
| 7 | §4.2 | KS test "D = 0.20, p > 0.85" vs Bateman–Horn | vs BH cumulative: D = 0.265, p = 0.62; vs uniform: D = 0.180, p = 0.95. | Minor |
| 8 | §5.2 | Poisson 1σ interval [3.4, 8.4]; P = 0.00781 | [3.4, 8.2] for mean 5.79 (now [3.5, 8.3] for 5.88); 8.94/(46×24.94) = 0.00779. | Trivial |
| 9 | §8.1 | "gap between p = 47 and p = 283 spans 61 primes" | 45 primes strictly between 47 and 283; 60 primes below 283 (all obstruction-free; 47 plays no special role). | Minor |
| 10 | §8.3 | "∼12,500 CPU-core-hours (∼2 weeks on 26 cores)" | 52 × 240 h = 12,500 h ignores the ≈1.3× cost of 600-digit arithmetic; ≈16,000 core-hours ≈ 4 weeks on 26 cores (12,500 h alone is already 20 days). | Minor |
| 11 | Eq. (1), all integrals | ln Q(t) ≈ 46 ln t | ln Q(t) = 46 ln t + ln 47 + O(1/t); the ln 47 term changes I₄ by −1.3 % and I₆ by −2.0 %. v2 keeps it. | Minor |
| 12 | §5 "no free parameters" | 𝔖₄ was fitted to the data (6,385) | v2 uses no calibration at all; the data-derived constant 749/I₄ = 6,533 agrees with the exact 6,514 to 0.3 %. | Framing |
| 13 | Repository | `paper/figures/` referenced by the .tex and README does not exist (figures are in `figures/`); `pdflatex` fails from `paper/`. `generate_figures.py` writes to `paper/figures/`. | v2 .tex uses `\graphicspath`; v2 figures in `figures/v2/`. | Repo |
| 14 | `compute_singular_series.py` | `milestones` are compared with `p in milestones`, so 100, 200, 500, … (non-primes) never print. | Fixed in v2 script (checkpoints flushed between primes). | Repo |

## Additions in v2

* **Cyclotomic identification** Q(n) = Φ₄₇(n, n−1) and the classical corollary that *every prime factor of every Q(n) is ≡ 1 (mod 47)* (Birkhoff–Vandiver). This gives a one-line proof of the bifurcated root structure and explains the 60-prime "unobstructed corridor".
* **Exact singular series** via 𝔖ₖ = [∏_{χ≠χ₀} L(1,χ)]^{−k} ∏ₚ Bₖ(p), with ∏ L(1,χ) = 0.1151693 (checked against L(1,χ₄₇) = π·5/√47) and an absolutely convergent product; values for k = 1…6 (`data/singular_series_exact.csv`).
* **Structure of ωₖ(p)**: ωₖ(p) = 46k for all but finitely many resonant primes; the exceptional primes below 10⁴ are listed (`data/resonant_primes_omega_v2.csv`, now k = 1…6).
* **Parameter-free hierarchy validation** (`data/hierarchy_validation.csv`): primes 18,473,571 vs 18,464,404 (+0.05 %), pairs 176,894 vs 176,778 (+0.07 %), triplets 1,779 vs 1,752, quadruplets 15 vs 17.0 (pioneer zone); 749 vs 746.9, 7 vs 5.88, 0 vs 0.047 (full range). Counting convention (maximal runs → positions) made explicit.
* **New Figure 5**: relative deviation of the plain partial product from the exact value for B ≤ 10⁶.
* Conditional model summed over the actual 742 quadruplets (E = 5.82) using `data/quadruplets_742.csv` copied from the Part II repository.
* Probability-of-detection column in the prediction table (50 % at ≈ 6.5×10¹², 63 % at N*, 90 % at ≈ 2.9×10¹³).

## Unchanged

The seven quintuplet starting values, the 742-quadruplet catalogue, Theorem 2.1 and its proof, the example p = 283 / n = 10, the ωₖ values at resonant primes, the qualitative emergence–extinction picture, and the order of magnitude of N*.

## Files

| v1 (kept unchanged) | v2 |
|---|---|
| `paper/Q47_Quintuplet_Boundary.tex/.pdf` | `paper/Q47_Quintuplet_Boundary_v2.tex/.pdf` |
| `scripts/compute_singular_series.py` | `scripts/compute_singular_series_v2.py` (exact + plain) |
| `scripts/generate_figures.py` | `scripts/generate_figures_v2.py`, `scripts/predictions_v2.py` |
| `figures/p2_fig*_v2.*` | `figures/v2/v2_fig1..5.*` |
| `data/singular_series_convergence.csv` (wrong rows) | `data/singular_series_convergence_v2.csv`, `data/singular_series_exact.csv`, `data/accelerated_convergence.csv`, `data/L1_product.txt` |
| `data/resonant_primes_omega.csv` | `data/resonant_primes_omega_v2.csv` (k = 1…6) |
| `data/quintuplets_7.csv` (digits −1, gap column shifted) | `data/quintuplets_7_v2.csv` |
| `data/sextuplet_predictions.csv`, `data/suppression_ratios.csv` | `data/sextuplet_predictions_v2.csv`, `data/suppression_ratios_v2.csv`, `data/hierarchy_validation.csv` |
| — | `data/quadruplets_742.csv` (copy of Part II) |

Reproduce everything with
```
python scripts/compute_singular_series_v2.py --bound 1000000 --write-csv
python scripts/predictions_v2.py
python scripts/generate_figures_v2.py
cd paper && pdflatex Q47_Quintuplet_Boundary_v2.tex && pdflatex Q47_Quintuplet_Boundary_v2.tex
```
