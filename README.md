# Q47 Quintuplet–Sextuplet Boundary

**Prime Quintuplets and the Sextuplet Boundary for Q(n) = n⁴⁷ − (n−1)⁴⁷**

*Part III of the Titan Project*

---

## Summary

Among the 742 prime quadruplets for Q(n) = n⁴⁷ − (n−1)⁴⁷ discovered over
1 ≤ n ≤ 2×10¹¹ ([Part II](https://github.com/Ruqing1963/Q47-Deep-Space-Quadruplet-Census)),
exactly **7 extend to prime quintuplets** (five consecutive n values all generating
probable primes of 487–519 digits), while **no sextuplet** was found.

| Constellation | Count | BH Prediction |
|:---|---:|---:|
| Quadruplets (k=4) | 742 | 742.0 (calibration) |
| **Quintuplets (k=5)** | **7** | **5.83** |
| Sextuplets (k=6) | 0 | 0.047 |

## Key Results

### Bifurcated Root Structure (Theorem 2.1)

Q(n) has a **bifurcated** root structure modulo primes:

- **ω₁(p) = 0** for all primes p ≢ 1 (mod 47) — no local obstruction
- **ω₁(p) = 46** for all primes p ≡ 1 (mod 47) — **resonant primes**

The resonant primes (p = 283, 659, 941, 1129, ...) impose a **periodic sieve
obstruction**, producing dramatic "cliffs" in the Bateman–Horn Euler product.
The first cliff at p = 283 halves the running singular series (σ₄(283) = 0.563).

### Corrected Singular Series

| Quantity | Value |
|:---|---:|
| 𝔖₄ (calibrated) | 6,385 |
| 𝔖₅ (corrected) | **57,108** |
| 𝔖₆ (corrected) | **519,756** |
| 𝔖₅/𝔖₄ | 8.94 |
| 𝔖₆/𝔖₅ | 9.10 |

### Sextuplet Boundary

The first sextuplet is predicted at:

> **N\* ≈ 1.05 × 10¹³**

requiring a 52-fold extension of the current search (~605-digit primes).

### Suppression Factor

The inter-level suppression factor is remarkably stable at **≈ 127×** across
five hierarchical levels spanning seven orders of magnitude—a structural
consequence of the Bateman–Horn framework.

## Repository Structure

```
├── paper/
│   ├── Q47_Quintuplet_Boundary.tex    # LaTeX source (10 pages)
│   ├── Q47_Quintuplet_Boundary.pdf    # Compiled paper
│   └── figures/
│       ├── p2_fig1_v2.{png,pdf}       # Quintuplet field + cumulative count
│       ├── p2_fig2_v2.{png,pdf}       # Morphological hierarchy + extension probs
│       ├── p2_fig3_v2.{png,pdf}       # Sextuplet prediction (log-log)
│       └── p2_fig4_v2.{png,pdf}       # Euler product cliff structure
├── data/
│   ├── quintuplets_7.csv              # All 7 quintuplet starting values
│   ├── resonant_primes_omega.csv      # ω_k(p) at all 28 resonant primes < 10000
│   ├── singular_series_convergence.csv # 𝔖_k(B) convergence table
│   ├── sextuplet_predictions.csv      # E[C₅], E[C₆] at various N
│   └── suppression_ratios.csv         # Inter-level suppression factors
├── scripts/
│   ├── compute_singular_series.py     # Reproduce 𝔖₄, 𝔖₅, 𝔖₆ computation
│   └── generate_figures.py            # Reproduce all 4 figures from data
├── README.md
└── LICENSE
```

## The Seven Quintuplets

| # | Starting n | Digits | n/10¹¹ |
|--:|---:|---:|---:|
| 1 | 35,676,017,721 | 487 | 0.357 |
| 2 | 64,482,563,907 | 498 | 0.645 |
| 3 | 73,292,417,435 | 501 | 0.733 |
| 4 | 116,255,850,744 | 510 | 1.163 |
| 5 | 147,743,683,226 | 515 | 1.477 |
| 6 | 159,430,471,996 | 516 | 1.594 |
| 7 | 182,501,065,420 | 519 | 1.825 |

## Reproducing Results

### Singular series computation
```bash
cd scripts
python compute_singular_series.py --bound 10000
```

### Regenerate figures
```bash
pip install matplotlib numpy scipy
cd scripts
python generate_figures.py
```

### Compile paper
```bash
cd paper
pdflatex Q47_Quintuplet_Boundary.tex
pdflatex Q47_Quintuplet_Boundary.tex
```

## Relation to Other Parts

| Part | Title | Repository |
|:---|:---|:---|
| I | Statistical Morphology (Pioneer Zone) | [Zenodo 18701355](https://zenodo.org/records/18701355) |
| II | Deep-Space Quadruplet Census | [GitHub](https://github.com/Ruqing1963/Q47-Deep-Space-Quadruplet-Census) / [Zenodo 18728540](https://zenodo.org/records/18728540) |
| **III** | **Quintuplets & Sextuplet Boundary** | **This repository** |

## Citation

```bibtex
@article{Chen2026c,
  author  = {Chen, Ruqing},
  title   = {Prime Quintuplets and the Sextuplet Boundary for
             {$Q(n) = n^{47} - (n-1)^{47}$}},
  year    = {2026},
  note    = {Part {III} of the Titan Project},
  url     = {https://github.com/Ruqing1963/Q47-Quintuplet-Sextuplet-Boundary}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Author

**Ruqing Chen** — GUT Geoservice Inc., Montréal, Canada
