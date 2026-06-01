from scipy.stats import norm
from math import ceil, log, sqrt

# ============================================================
# Complexity Estimator for the Sparsity Below Threshold Algorithm
# ============================================================
#
# This script implements the complexity analysis presented
# in Section 4 for ternary sLWE instances whose sparsity
# parameter lies below the threshold limit.
#
# For each parameter set (n,k,m), the algorithm:
#
#   1. Estimates the size of a subsystem covering all
#      variables.
#
#   2. Estimates the cost of constructing such a subsystem.
#
#   3. Estimates the cost of solving the subsystem.
#
#   4. Reports the resulting time and memory exponents.
#
# All complexity values are reported as base-2 exponents.
#
# ============================================================


def complexity_estimator(n, k, m, alpha):

    # ========================================================
    #  Subsystem Construction
    # ========================================================
    #
    # Let V_i denote the number of distinct variables covered
    # after selecting i equations.
    #
    # Initially, a single equation contributes k variables.
    # The subsystem-construction procedure iteratively adds
    # equations until all variables are covered.
    #
    # The coverage process is approximated using the mean and
    # variance derived in Theorem 1.
    #
    # ========================================================

    covered_variables = k
    subsystem_size = 0

    while covered_variables < n:

        subsystem_size += 1

        # ----------------------------------------------------
        # Coverage Distribution
        # ----------------------------------------------------
        #
        # Let X denote the number of variables covered after
        # adding one additional equation to a subsystem
        # currently covering covered_variables variables.
        #
        # The following expressions compute
        #
        #     E[X]
        #
        # and
        #
        #     Var[X].
        #
        # ----------------------------------------------------

        mean = (
            n
            * (
                1
                - (1 - covered_variables / n)
                * (1 - k / n)
            )
        )

        variance = (
            n
            * (1 - covered_variables / n)
            * (1 - k / n)
            - n**2
            * (1 - covered_variables / n) ** 2
            * (1 - k / n) ** 2
            + n
            * (n - 1)
            * (
                1
                - 2 * covered_variables / n
                + covered_variables
                * (covered_variables - 1)
                / (n * (n - 1))
            )
            * (
                1
                - 2 * k / n
                + k * (k - 1)
                / (n * (n - 1))
            )
        )

        sd = sqrt(variance)

        # ----------------------------------------------------
        # Conservative Coverage Estimate
        # ----------------------------------------------------
        #
        # The next subsystem size is estimated by
        #
        #     mean + alpha * sd.
        #
        # This corresponds to moving alpha standard
        # deviations above the mean.
        #
        # ----------------------------------------------------

        next_coverage = ceil(mean + alpha * sd)

        # ----------------------------------------------------
        # Tail Probability Estimate
        # ----------------------------------------------------
        #
        # The normal approximation is used to estimate the
        # probability that the number of covered variables
        # exceeds the chosen threshold.
        #
        # The reciprocal of this probability provides a
        # heuristic estimate of the number of trials
        # required to obtain such an equation.
        #
        # ----------------------------------------------------

        tail_probability = (
            1
            - norm.cdf(
                next_coverage - 0.5,
                loc=mean,
                scale=sd
            )
        )

        expected_trials = 1 / tail_probability

        # ----------------------------------------------------
        # Feasibility Condition
        # ----------------------------------------------------
        #
        # If the expected number of trials exceeds the total
        # number of available equations, the subsystem
        # construction process is considered infeasible.
        #
        # ----------------------------------------------------

        if m < expected_trials:

            subsystem_size += (n - covered_variables)

            print(
                f"(n,k)=({n},{k}) "
                f"Subsystem size ≈ {subsystem_size}"
            )

            break

        if next_coverage == covered_variables:

            print(
                f"(n,k)=({n},{k}) "
                f"Subsystem size ≈ {subsystem_size}"
            )

            break

        covered_variables = next_coverage

    # ========================================================
    # Solving the Constructed Subsystem
    # ========================================================
    #
    # Once a subsystem covering all variables has been
    # obtained, the complexity of solving the subsystem is
    # estimated using the analysis of Section 4.
    #
    # Let
    #
    #     kappa = ceil((k+1)/2).
    #
    # ========================================================

    kappa = ceil((k + 1) / 2)

    # --------------------------------------------------------
    # Time Complexity Components
    # --------------------------------------------------------
    #
    # T_construct:
    #     Cost of constructing the subsystem.
    #
    # T_solve:
    #     Cost of solving the subsystem equations.
    #
    # --------------------------------------------------------

    T_construct = subsystem_size * m

    T_solve = subsystem_size * (3 ** kappa)

    # ========================================================
    # Overall Time and Memory Complexity
    # ========================================================
    #
    # Time:
    #
    #     max{
    #         T_construct,
    #         T_solve
    #     }.
    #
    # Memory:
    #
    #     max{
    #         3^kappa,
    #         m
    #     }.
    #
    # All values are reported as base-2 exponents.
    #
    # ========================================================

    time_exponent = ceil(
        log(max(T_construct, T_solve), 2)
    )

    memory_exponent = ceil(
        log(max(3 ** kappa, m), 2)
    )

    # ========================================================
    # Output
    # ========================================================

    print("=" * 70)
    print("PARAMETERS")
    print("=" * 70)

    print(f"n      = {n}")
    print(f"k      = {k}")
    print(f"m      = {m}")
    print(f"alpha  = {alpha}")

    print()

    print("=" * 70)
    print("COMPLEXITY EXPONENTS")
    print("=" * 70)

    print(f"Subsystem size estimate     = {subsystem_size}")
    print(f"log2(Time)                 = {time_exponent}")
    print(f"log2(Memory)               = {memory_exponent}")

    print()
    print("=" * 70)
    print()


# ============================================================
# Parameter Sets from Table 1 and Table 4
# ============================================================
#
# Each tuple has the form
#
#     (n, m, k, alpha)
#
# and corresponds to one row of the complexity table.
#
# ============================================================

params = [

    (1218, 2**13, 20, 2),
    (1143, 2**13, 30, 3),
    (1110, 2**13, 40, 3),

    (1425, 2**17, 20, 3),
    (1265, 2**17, 30, 3),
    (1195, 2**17, 40, 3),

    (1656, 2**21, 20, 4),
    (1395, 2**21, 30, 4),
    (1285, 2**21, 40, 4),
]

for (n, m, k, alpha) in params:
    complexity_estimator(n, k, m, alpha)