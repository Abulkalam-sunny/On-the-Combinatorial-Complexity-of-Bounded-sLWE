This repository contains the source code accompanying the paper:

**“On the Combinatorial Complexity of Bounded Sparse Learning with Errors (sLWE)”**

The repository provides implementations of the complexity estimators, experimental verification procedures, and numerical validation scripts used to generate the results reported throughout the paper.

---

## Repository Structure

```
├── Complexity_Estimators/
│   ├── complexity_estimator_below_threshold.py
│   └── complexity_estimator_above_threshold.py
│
├── Experimental_Verification/
│   ├── experimental_verification_algorithm4.sage
│   └── experimental_verification_algorithm5.sage
│
└── Numerical_Validation/
    └── normal_approximation_validation.py
    
```

---

## Description of Files

### Complexity_Estimators/

#### `complexity_estimator_below_threshold.py`

Implements the heuristic complexity analysis for solving ternary sLWE instances whose sparsity parameter lies **below the threshold limit**.

This script estimates:

* subsystem construction complexity,
* subsystem solving complexity,
* overall running time,
* memory complexity.

The script reproduces the complexity exponents reported in Table~4 of the paper.

---

#### `complexity_estimator_above_threshold.py`

Implements the heuristic complexity analysis for solving ternary sLWE instances whose sparsity parameter lies **above the threshold limit**.

This implementation includes:

* subsystem extraction,
* subsystem solving,
* subsystem merging,
* overall time estimation,
* memory estimation,
* feasibility checks.

The script reproduces the complexity estimates reported in Table 5.

---

### Experimental_Verification/

#### `experimental_verification_algorithm4.sage`

Experimental verification of **Algorithm 4**.

The implementation generates random ternary sLWE instances and empirically evaluates:

* subsystem construction,
* variable coverage,
* intermediate table sizes,
* recovery performance.

---

#### `experimental_verification_algorithm5.sage`

Experimental verification of **Algorithm 5**.

The implementation experimentally evaluates:

* subsystem extraction,
* iterative subsystem solving,
* subsystem merging,
* final solution recovery.

The script reproduces the experimental statistics reported in Section 7.

---

### Numerical_Validation/

#### `normal_approximation_validation.py`

Provides numerical validation of the normal approximation used in subsystem extraction.

The script:

* computes theoretical means and variances,
* performs Monte Carlo simulation,
* validates concentration within 1σ, 2σ, and 3σ intervals,
* produces the distribution figures reported in Appendix D.

---

## Requirements

The scripts were tested using:

```text
Python ≥ 3.10
SageMath ≥ 10
NumPy
SciPy
Matplotlib
```

Install Python dependencies using:

```bash
pip install numpy scipy matplotlib
```

---

## Reproducibility

All scripts are self-contained and reproduce the complexity estimates, experimental statistics, and numerical validation reported in the paper.

Randomized experiments use independently generated ternary sLWE instances and therefore minor numerical variation may occur across executions.

---

