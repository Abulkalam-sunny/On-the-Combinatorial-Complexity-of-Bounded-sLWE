# ============================================================
# Experimental Verification of the Algorithm for Solving
# Ternary sLWE Instances with Sparsity Below the Threshold
# Limit
# ============================================================


# This script implements the proof-of-concept experiments
# described in Section 4.


# Experimental workflow:


# Phase I : Generate random ternary sLWE instances
# Phase II : Construct a covering subsystem
# Phase III : Rearrange equations for maximal overlap
# Phase IV : Solve and merge equations iteratively
# Phase V : Evaluate uniqueness statistics


# Statistics are averaged over repeated experiments.


# ============================================================



import itertools
from random import sample
from collections import defaultdict
import random

# ==============================================================
# Solve Algorithm : Section : 3.2
# ============================================================
    # Solve a sparse ternary equation using a meet-in-the-middle
    # enumeration strategy.

    # Given

    #     <y,s> = b (mod M),

    # where s ∈ {-1,0,1}^n, this function enumerates all ternary
    # assignments satisfying the equation.

    # Parameters
    # ----------
    # y : vector
    #     Sparse coefficient vector.

    # b : integer
    #     Right-hand-side value modulo M.

    # Returns
    # -------
    # list
    #     [solution_table, active_variables]

    #     solution_table :
    #         List of ternary solution vectors.

    #     active_variables :
    #         Indices of variables appearing in y.

# ============================================================


def solve(y, b):
    # ========================================================
    # Extract Active Variables
    # ========================================================
    #
    # Remove zero coefficients and record only the variables
    # participating in the equation.
    #
    # ========================================================

    active_variables = []
    active_coefficients = []

    for position in range(len(y)):

        if y[position] != 0:

            active_variables.append(position)
            active_coefficients.append(y[position])

    num_variables = len(active_coefficients)

    # ========================================================
    # Split Variables (Meet-in-the-Middle)
    # ========================================================
    #
    # Divide active variables into two halves.
    #
    # ========================================================

    left_coefficients = active_coefficients[: num_variables // 2]

    right_coefficients = active_coefficients[num_variables // 2 :]

    # ========================================================
    # Enumerate Left Half
    # ========================================================
    #
    # Store all partial sums in a hash table indexed by
    # modular inner products.
    #
    # ========================================================

    left_table = defaultdict(list)

    for left_assignment in itertools.product([-1, 0, 1],repeat=len(left_coefficients)):

        partial_sum = vector(Zmod(M), list(left_assignment)) * vector(Zmod(M), left_coefficients)

        left_table[int(partial_sum)].append(left_assignment)

    # ========================================================
    # Enumerate Right Half and Merge
    # ========================================================

    solution_table = []

    for right_assignment in itertools.product([-1, 0, 1], repeat=len(right_coefficients)):

        remaining_value = b - vector(Zmod(M), list(right_assignment)) * vector(Zmod(M), right_coefficients)
        

        for err in [-1, 0, 1]:

            lookup_value = int( remaining_value + err)

            if left_table[lookup_value]:

                for left_assignment in left_table[lookup_value]:
                    
                    full_solution = [0] * len(y)

                    merged_assignment =  list(left_assignment) + list(right_assignment)

                    for position in range(len(active_variables)):

                        full_solution[active_variables[position]] = merged_assignment[position]
                        

                    solution_table.append(vector(Zmod(M), full_solution))

    # ========================================================
    # Output
    # ========================================================

    return [solution_table, active_variables]



# ==============================================================
# Merge Algorithm : Section : 3.3
# ============================================================
    # Merge two solution tables by matching assignments on
    # common variables.

    # Let

    #     S1 = (T1, V1)
    #     S2 = (T2, V2)

    # denote two solution tables together with their active
    # variable sets.

    # Solutions are merged by requiring consistency on the
    # variables appearing in both subsystems.

    # Parameters
    # ----------
    # left_table : list
    #     [solutions_1, variables_1]

    # right_table : list
    #     [solutions_2, variables_2]

    # Returns
    # -------
    # list
    #     [merged_solution_table, merged_variables]

    #     merged_solution_table :
    #         Combined compatible assignments.

    #     merged_variables :
    #         Union of active variables.

# ============================================================

def Merge(left_table, right_table):
    # ========================================================
    # Extract Tables and Variable Sets
    # ========================================================

    solutions_left, variables_left = left_table
    solutions_right, variables_right = right_table

    # ========================================================
    # Handle Degenerate Cases
    # ========================================================
    #
    # If one subsystem contains no active variables, return
    # the remaining subsystem unchanged.
    #
    # ========================================================

    if len(variables_left) == 0:
        return [solutions_right, variables_right]

    if len(variables_right) == 0:
        return [solutions_left,variables_left]

    # ========================================================
    # Compute Common and Merged Variable Sets
    # ========================================================
    #
    # Common variables determine compatibility conditions.
    #
    # ========================================================

    common_variables = list(set(variables_left).intersection(set(variables_right)))

    merged_variables = list(set(variables_left).union(set(variables_right)))

    # ========================================================
    # Hash Left Table by Projection
    # ========================================================
    #
    # Each solution is indexed by its restriction to the
    # common variables.
    #
    # ========================================================

    projection_table = {}

    for solution_left in solutions_left:

        projection = tuple(solution_left[position] for position in common_variables)

        projection_table.setdefault(projection,[]).append(solution_left)

    # ========================================================
    # Probe Right Table and Merge
    # ========================================================
    #
    # Solutions are compatible if their assignments agree on
    # all common variables.
    #
    # ========================================================

    merged_solution_table = []

    for solution_right in solutions_right:

        projection = tuple(solution_right[position] for position in common_variables)

        if projection not in projection_table:
            continue

        for solution_left in projection_table[projection]:

            merged_assignment = [0] * n

            # --------------------------------------------
            # Copy coordinates from left solution
            # --------------------------------------------

            for position in variables_left:

                merged_assignment[position] = solution_left[position]

            # --------------------------------------------
            # Copy coordinates from right solution
            # --------------------------------------------

            for position in variables_right:

                merged_assignment[position] = solution_right[position]

            merged_solution_table.append(vector(Zmod(M),merged_assignment))

    # ========================================================
    # Output
    # ========================================================

    return [merged_solution_table,merged_variables]







# ============================================================
# Experimental Verification:
# Ternary sLWE with Sparsity Below the Threshold Limit
# ============================================================

subsystem_size_statistics = 0

uniqueness_statistics = 0

uniqueness_distribution = []

num_experiments = 1000


# ============================================================
# Main Experimental Loop
# ============================================================

for experiment in range(1,num_experiments + 1):

    import random

    # ========================================================
    # Generate Random Ternary sLWE Instance
    # ========================================================

    n = 1218
    m = 2**13

    M = 2**32

    k = 20

    coefficient_matrix = Matrix(Zmod(M),m,n,lambda i, j: 0)

    support_sets = []

    for row in range(m):

        support = random.sample(range(n),k)

        support_sets.append(support)

        for position in support:

            coefficient_matrix[row,position] = randint(1,M - 1)

    secret = vector(Zmod(M),[random.choice([-1, 0, 1]) for _ in range(n)])

    error = vector(Zmod(M),[random.choice([-1, 0, 1]) for _ in range(m)])

    target = coefficient_matrix * secret - error

    # ========================================================
    # Construct Covering Subsystem
    # ========================================================

    initial_rows = random.sample(range(m),1)

    covered_variables = set()

    for row in initial_rows:

        covered_variables.update(support_sets[row])

    subsystem_rows = list(initial_rows)

    while len(covered_variables) < n:

        best_gain = -1

        for row in range(m):

            gain = len(set(support_sets[row])-covered_variables)

            if gain > best_gain:

                best_gain = gain

                best_row = row

        covered_variables.update(support_sets[best_row])

        subsystem_rows.append(best_row)

    subsystem_size_statistics += len(subsystem_rows)

    # ========================================================
    # Overlap-Based Equation Reordering
    # ========================================================

    reordered_rows = [subsystem_rows[0]]

    overlap_variables = set(support_sets[subsystem_rows[0]])

    while (len(reordered_rows)<len(subsystem_rows)):

        best_overlap = -1

        for row in (set(subsystem_rows) - set(reordered_rows)):

            overlap = len(set(support_sets[row]) & overlap_variables)

            if overlap > best_overlap:

                best_overlap = overlap

                best_row = row

        reordered_rows.append(best_row)

        overlap_variables.update(support_sets[best_row])

    # ========================================================
    # Solve and Merge Equations
    # ========================================================

    merged_table, merged_variables = solve(coefficient_matrix[reordered_rows[0]],target[reordered_rows[0]])

    

    uniqueness_equation_count = 1

    equations_used = 1

    for row in reordered_rows[1:]:

        equations_used += 1

        current_table, current_variables = solve(coefficient_matrix[row],target[row])



        merged_table, merged_variables = Merge((merged_table,merged_variables),(current_table,current_variables))



        # ----------------------------------------------------
        # Record the last equation for which multiple
        # candidate solutions remain.
        # ----------------------------------------------------

        if len(merged_table) > 1:

            uniqueness_equation_count = equations_used + 1


    # ========================================================
    # Aggregate Statistics
    # ========================================================

    uniqueness_distribution.append(uniqueness_equation_count)

    uniqueness_statistics += (uniqueness_equation_count)

    print("=" * 60)

    print("Experiment",experiment,"| Success:",secret in merged_table,'\n')


    print("Average equations for coverage =",round(subsystem_size_statistics/experiment,2))

    print("Average equations for uniqueness =",round(uniqueness_statistics/experiment,2))


    print()

# ============================================================
# Final Statistics
# ============================================================

print("_" * 60)

print("Maximum equations required for uniqueness =",max(uniqueness_distribution))
