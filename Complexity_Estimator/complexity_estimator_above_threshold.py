import numpy as np
import math
from math import comb
from scipy.stats import norm

# ============================================================
# Auxiliary Mathematical Functions
# ============================================================
#
# This contains helper routines used throughout the
# complexity estimator.
#
# log2_binomial(n,k):
#     Computes log2(binomial(n,k)) using the logarithmic gamma
#     function. This avoids overflow and allows feasibility
#     checks to be performed for large parameter sets.
#
# ============================================================

def log2_binomial(n,k):
    return (
        math.lgamma(n+1)
        - math.lgamma(k+1)
        - math.lgamma(n-k+1)
    ) / math.log(2)

# ============================================================
# Dynamic Programming Procedure from Appendix B
# ============================================================
#
# Let X_t denote the number of distinct variables appearing
# after selecting t number of k-sparse equations.
#
# The function extradraws() computes the probability
# distribution of X_t using the procedure
# described in Appendix B.
#
# This routine is subsequently used to evaluate conditional
# expectations of the form
#
#     E[X_gamma | X_mu = v]
#
# which arise in the analysis of subsystem solving.
#
# ============================================================

def extradraws(numberextradraws,
               initial,
               total,
               draw):

    probunique = np.zeros(total)
    probunique[initial-1] = 1.0

    for _ in range(numberextradraws):

        newprobunique = np.zeros(total)

        for already in range(total):

            if probunique[already] == 0:
                continue

            for future in range(already,total):

                overlap = draw-(future-already)

                if overlap < 0:
                    continue

                num1 = comb(
                    total-already-1,
                    future-already
                )

                num2 = comb(
                    already+1,
                    overlap
                )

                den = comb(total,draw)

                prob = num1*num2/den

                newprobunique[future] += (
                    probunique[already]*prob
                )

        probunique = newprobunique

    return probunique

# ============================================================
# Complexity Estimator for the Sparsity Above Threshold Algorithm
# ============================================================
#
# Input:
#     n  : number of variables
#     k  : sparsity parameter
#     m  : number of equations
#     p  : modulus exponent (M = 2^p)
#     mu : subsystem size
#     v  : subsystem variable bound
#     r  : number of subsystems
#
# Output:
#     Estimates of:
#
#         T_get_subsystem
#         T_solve_subsystem
#         T_merge
#         Total running time
#         Memory complexity
#
# All complexity values are reported as base-2 exponents.
#
# ============================================================


def complexity_estimator(n,k,m,p,mu,v,r):
    delta = int(p/math.log(3,2))


    # ============================================================
    # Subsystem Extraction Cost
    # ============================================================
    #
    # The number of variables appearing in a random collection of
    # mu equations is approximated by a normal distribution whose
    # mean and variance are derived in Theorem 1.
    #
    # The probability that a randomly selected subsystem contains
    # at most v variables is used to estimate the expected number
    # of trials required to obtain a suitable subsystem.
    #
    #     T_get_subsystem ≈ 1 / Pr(X_mu ≤ v).
    #
    # ============================================================
    
    mean = n*(1-(1-k/n)**mu)
    
    var = (
        n*(1-k/n)**mu
        - n**2*(1-k/n)**(2*mu)
        + n*(n-1)
          * (
              1
              - 2*k/n
              + k*(k-1)/(n*(n-1))
            )**mu
    )
    
    sd = math.sqrt(var)
    
    prob = norm.cdf(
        v + 0.5,
        loc = mean,
        scale = sd
    )
    
    E_get = -math.log2(prob)
    
    # ============================================================
    # Feasibility Condition
    # ============================================================
    #
    # The expected number of subsystem searches required to
    # obtain r suitable subsystems must not exceed the total
    # number of available choices of mu equations.
    #
    # The condition checked here is
    #
    #     r · T_get_subsystem < binomial(m, mu).
    #
    # ============================================================
    
    E_get_total = math.log2(r) + E_get
    
    feasible = (
        E_get_total
        <
        log2_binomial(m,mu)
    )
    
    # ============================================================
    # Cost of Solving a Single Equation
    # ============================================================
    #
    # Implements the complexity estimate derived in Section 3.
    #
    # The cost is taken as the larger of the meet-in-the-middle
    # estimate and the collision-based estimate, expressed as a
    # base-2 exponent.
    #
    # ============================================================
    
    single_equation_exp = max(
        math.ceil((k+1)/2)*math.log2(3),
        (k+1)*math.log2(3)-p
    )
    
    # ============================================================
    # Solving a Subsystem
    # ============================================================
    #
    # For each gamma = 1,...,mu, the algorithm computes
    #
    #     E[X_gamma | X_mu = v]
    #
    # using the procedure from Appendix B.
    #
    # These conditional expectations are then used to estimate
    # the sizes of intermediate solution tables occurring during
    # iterative equation merging inside a subsystem.
    #
    # The largest intermediate table determines the estimate
    #
    #     T_solve_subsystem.
    #
    # ============================================================
    
    max_intermediate_exp = -10**100
    
    for gam in range(1,mu+1):
    
        probs_after_gam = extradraws(
            numberextradraws = gam-1,
            initial = k,
            total = n,
            draw = k
        )
    
        numerator = np.zeros(n)
    
        for covered_after_gam in np.where(
            probs_after_gam > 0
        )[0]:
    
            tail = extradraws(
                numberextradraws = mu-gam,
                initial = covered_after_gam+1,
                total = n,
                draw = k
            )
    
            numerator[covered_after_gam] = (
                probs_after_gam[covered_after_gam]
                *
                tail[v-1]
            )
    
        if np.sum(numerator) == 0:
            continue
    
        conditional_distribution = (
            numerator /
            np.sum(numerator)
        )
    
        expected_value_given_v = np.sum(
            np.arange(1,n+1)
            *
            conditional_distribution
        )
    
        exponent = (
            math.ceil(
                (
                    expected_value_given_v
                    + gam
                )
                *
                math.log2(3)
            )
            -
            p*gam
        )
    
        max_intermediate_exp = max(
            max_intermediate_exp,
            exponent
        )
    
    E_solve_subsystem = max(
        math.log2(mu) + single_equation_exp,
        max_intermediate_exp
    )
    
    E_solve_total = (
        math.log2(r)
        +
        E_solve_subsystem
    )
    
    # ============================================================
    # Merging Subsystems
    # ============================================================
    #
    # Let v_j denote the expected number of variables appearing
    # after merging j subsystems.
    #
    # The expected merged-table size is estimated by
    #
    #     3^(v_j + j*mu) / M^(j*mu).
    #
    # The largest such quantity over all
    #
    #     j = 1,...,r
    #
    # determines T_merge.
    #
    # ============================================================
    
    E_merge = -10**100
    
    for j in range(1,r+1):
    
        vj = n * (
            1
            -
            (1-v/n)**j
        )
    
        exponent = (
            (vj + j*mu)
            *
            math.log2(3)
            -
            p*j*mu
        )
    
        E_merge = max(
            E_merge,
            exponent
        )
    
    # ============================================================
    # Overall Time Complexity
    # ============================================================
    #
    # The total running time is estimated as
    #
    #     max{
    #         r*T_get_subsystem,
    #         r*T_solve_subsystem,
    #         T_merge
    #     }.
    #
    # All quantities are reported as base-2 exponents.
    #
    # ============================================================
    
    E_total = max(
        E_get_total,
        E_solve_total,
        E_merge
    )
    
    # ============================================================
    # OUTPUT
    # ============================================================
    
    print("="*70)
    print("PARAMETERS")
    print("="*70)
    
    print(f"n      = {n}")
    print(f"k      = {k}")
    print(f"m      = {m}")
    print(f"p      = {p}")
    print(f"mu     = {mu}")
    print(f"v      = {v}")
    print(f"r      = {r}")
    
    print()
    print("="*70)
    print("EXPONENTS")
    print("="*70)
    
    print(f"log2(T_get_subsystem)       = {E_get:.4f}")
    print(f"log2(r*T_get_subsystem)     = {E_get_total:.4f}")
    
    print()
    
    print(f"log2(T_solve_subsystem)     = {E_solve_subsystem:.4f}")
    print(f"log2(r*T_solve_subsystem)   = {E_solve_total:.4f}")
    
    print()
    
    print(f"log2(T_merge)               = {E_merge:.4f}")
    
    print()
    
    print(f"log2(T_total)               = {E_total:.4f}")

    # ============================================================
    # Memory Complexity
    # ============================================================
    #
    # The memory requirement is estimated as the largest of:
    #
    #     T_merge,
    #     T_solve_subsystem,
    #     m.
    #
    # This corresponds to the largest table that must be stored
    # during the attack.
    #
    # ============================================================
    
    print(f"log2(Memory)               = {max(E_merge,E_solve_subsystem,math.log(m,2)):.4f}")
    
    print()
    print("="*70)
    print("FEASIBILITY")
    print("="*70)
    
    print(
        f"log2(r*T_get_subsystem) = {E_get_total:.4f}"
    )
    
    print(
        f"log2(binomial(m,mu)) = {log2_binomial(m,mu):.4f}"
    )
    
    print(
        f"Feasible = {feasible}"
    )
    print("="*70)
    print("="*70)
    print('\n')


# ============================================================
# Parameter Sets from Table 5
# ============================================================
#
# Each tuple has the form
#
#     (n, k, m, p, mu, v, r)
#
# and corresponds to one row of Table 5.
#
# ============================================================

params = [

    # =====================================================
    # M = 2^64
    # =====================================================

    (1090, 50, 2**13, 64,  2,  92, 8),
    (1158, 50, 2**17, 64,  2,  92, 8),
    (1234, 50, 2**21, 64,  2,  92, 9),

    # =====================================================
    # M = 2^32
    # =====================================================

    (1143, 30, 2**13, 32,  8, 186, 4),
    (1265, 30, 2**17, 32,  9, 208, 4),
    (1395, 30, 2**21, 32,  9, 209, 4),

    (1110, 40, 2**13, 32, 11, 289, 4),
    (1195, 40, 2**17, 32, 11, 288, 4),
    (1285, 40, 2**21, 32, 12, 313, 4),

    (1090, 50, 2**13, 32, 13, 387, 4),
    (1158, 50, 2**17, 32, 13, 375, 4),
    (1234, 50, 2**21, 32, 13, 373, 4),

    # =====================================================
    # M = 2^16
    # =====================================================

    (1218, 20, 2**13, 16, 24, 302, 4),
    (1425, 20, 2**17, 16, 28, 352, 4),
    (1656, 20, 2**21, 16, 33, 414, 4),

    (1143, 30, 2**13, 16, 28, 421, 4),
    (1265, 30, 2**17, 16, 30, 448, 4),
    (1395, 30, 2**21, 16, 33, 492, 4),

    (1110, 40, 2**13, 16, 29, 528, 4),
    (1195, 40, 2**17, 16, 30, 515, 4),
    (1285, 40, 2**21, 16, 32, 546, 4),

    (1090, 50, 2**13, 16, 29, 611, 4),
    (1158, 50, 2**17, 16, 30, 600, 5),
    (1234, 50, 2**21, 16, 31, 591, 5),
]

for (n,k,m,p,mu,v,r) in params:
    complexity_estimator(n,k,m,p,mu,v,r)