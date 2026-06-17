# bezier_shuttle
Reproducibility code accompanying the manuscript:
"A Bézier Shuttle for Accelerating Brent’s Root-Finder" 

## Quick start

```bash
git clone https://github.com/oswcad/bezier_shuttle.git
cd bezier_shuttle
pip install -r requirements.txt
python table1.py
```

# Bézier Shuttle for Brent's Root Finder

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()

Reproducibility code for the paper:

> **A Bézier Shuttle for Accelerating Brent’s Root-Finder**  
> Oswaldo Cadenas  


This repository reproduces **Table 1** of the paper, comparing the original 1973 Brent method with the proposed two-step Bézier Shuttle + Brent hybrid.  
The experiments confirm that the Bézier Shuttle safely reduces the iteration count for curved and oscillatory functions while preserving bracketing stability and final residual accuracy.

---

## 📦 Requirements

- Python ≥ 3.9  
- NumPy ≥ 1.23

You can install the dependencies with:

```bash
pip install -r requirements.txt

```
## Usage
To reproduce Table 1:
```bash
python table1.py
```

## Expected output
```
Function                             α      Brent    Hybrid   Gain        Root
--------------------------------------------------------------------------------
e^(-x+11)-2                          0.5       40       30     10    10.306853
0.1(x-2)(x-5)(x-8)                   0.5        1        1      0     5.000000
sin(x)+0.5x-2                        0.5       27       27      0     5.462807
(x-3.5)^(-1)+2                       0.5       27       21      6     3.000000
ln(x+1)-2                            0.5       39       28     11     6.389056
x^3-5x+1                        adaptive       25       13     12     2.128419
cos(x)-x                             0.5       11       13     -2     0.739085
e^(-x)-0.1                           0.5       37       35      2     2.302585
```

## Background
	•	Brent (1973): Algorithms for Minimization Without Derivatives — the classical hybrid bisection/secant/IQI method.
	•	Bézier Shuttle (Cadenas 2025): a two-step quadratic Bézier contraction that precedes Brent’s landing phase. Each shuttle uses the midpoint ( x_c ) to curve the bracket toward the function’s shape before invoking Brent’s interpolation.

## Citation

If you use this software, please cite the accompanying manuscript:

> Oswaldo Cadenas,

> *A Bézier Shuttle for Accelerating Brent's Root-Finder*,

> manuscript under review.

## License
MIT License © 2025 Oswaldo Cadenas

---

## 🧰 2. requirements.txt

```txt
numpy>=1.23

