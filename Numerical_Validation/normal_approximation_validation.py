# ============================================================
# Numerical Validation of the Normal Approximation
# for Subsystem Extraction
# ============================================================
#
# This script provides numerical evidence supporting the
# normal approximation used in Appendix D.
#
# Let X_mu denote the number of distinct variables appearing
# after selecting mu sparse equations uniformly at random.
#
# For each parameter set:
#
#   1. Compute the theoretical mean and variance.
#   2. Generate Monte Carlo samples.
#   3. Compare empirical and theoretical statistics.
#   4. Validate 1σ, 2σ, and 3σ concentration behavior.
#   5. Plot empirical distributions against the normal
#      approximation.
#
# ============================================================


# ============================================================
# Experimental Parameters
# ============================================================


import random
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# ============================================================
# Parameters
# ============================================================

m = 2**17
Nexp = 100000

param_sets = [
    (1425, 20),
    (1265, 30),
    (1195, 40),
    (1158, 50)
]

mu_values = [10, 20, 30, 40]

# ============================================================
# Plot
# ============================================================

fig, axes = plt.subplots(
    len(param_sets),
    len(mu_values),
    figsize=(20, 18)
)

for row, (n, k) in enumerate(param_sets):

    print("\n" + "=" * 100)
    print(f"(n,k)=({n},{k})")
    print("=" * 100)

    # Generate sparse supports once for this (n,k)
    Mat_idx = [
        random.sample(range(n), k)
        for _ in range(m)
    ]

    for col, mu in enumerate(mu_values):

        # ----------------------------------------------------
        # Theoretical mean
        # ----------------------------------------------------

        mean = n * (1 - (1 - k/n)**mu)

        # ----------------------------------------------------
        # Theoretical variance
        # ----------------------------------------------------

        var = (
            n * (1 - k/n)**mu
            - n**2 * (1 - k/n)**(2 * mu)
            + n * (n - 1)
            * (
                1
                - 2 * k / n
                + k * (k - 1) / (n * (n - 1))
            ) ** mu
        )

        sd = math.sqrt(var)

        # ----------------------------------------------------
        # Monte Carlo experiment
        # ----------------------------------------------------

        samples = []

        for _ in range(Nexp):

            idx = random.sample(range(m), mu)

            S = set()

            for i in idx:
                S.update(Mat_idx[i])

            samples.append(len(S))

        samples = np.array(samples)

        # ----------------------------------------------------
        # Empirical statistics
        # ----------------------------------------------------

        emp_mean = samples.mean()
        emp_var = samples.var()

        # ----------------------------------------------------
        # Normal approximation validation
        # ----------------------------------------------------

        emp1 = np.mean(np.abs(samples - mean) <= sd)
        emp2 = np.mean(np.abs(samples - mean) <= 2 * sd)
        emp3 = np.mean(np.abs(samples - mean) <= 3 * sd)

        th1 = norm.cdf(1) - norm.cdf(-1)
        th2 = norm.cdf(2) - norm.cdf(-2)
        th3 = norm.cdf(3) - norm.cdf(-3)

        print(
            f"mu={mu:2d} | "
            f"Mean: {emp_mean:.2f}/{mean:.2f} | "
            f"Var: {emp_var:.2f}/{var:.2f} | "
            f"1σ: {emp1:.4f}/{th1:.4f} | "
            f"2σ: {emp2:.4f}/{th2:.4f} | "
            f"3σ: {emp3:.4f}/{th3:.4f}"
        )

        # ----------------------------------------------------
        # Plot
        # ----------------------------------------------------

        ax = axes[row, col]

        # Histogram
        ax.hist(
            samples,
            bins=range(min(samples), max(samples)+2),
            density=True,
            alpha=0.7
        )

        # Fitted normal curve
        x = np.linspace(
            samples.min() - 2,
            samples.max() + 2,
            1000
        )

        ax.plot(
            x,
            norm.pdf(x, mean, sd),
            linewidth=2
        )

        ax.set_title(
            rf"$n={n},\,k={k},\,\mu={mu}$",
            fontsize=11
        )

        ax.set_xlabel("Number of variables")
        ax.set_ylabel("Density")

# ============================================================
# Layout
# ============================================================

plt.tight_layout()
plt.savefig("normal_validation.pdf", dpi=300)