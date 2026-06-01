# Experimental Verification of the Algorithm for Solving
# Ternary sLWE Instances with Sparsity Above the Threshold
# Limit
# ============================================================


# This script implements the proof-of-concept experiments
# described in Section 5.


# Experimental workflow:


# Phase I : Generate random ternary sLWE instances
# Phase II : Extract sparse subsystems
# Phase III : Rearrange equations for maximal overlap
# Phase IV : Solve individual subsystems
# Phase V : Merge subsystem solution tables
# Phase VI : Recover remaining variables
# Phase VII : Aggregate statistics


# Reported quantities:


# • subsystem extraction effort
# • intermediate variable counts
# • intermediate solution-table sizes
# • additional equations required


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
# Experimental Parameters
# ============================================================

n = 220
m = 2**10

M = 2**16

k = 16

v = 47
mu = 4
r = 4

num_experiments = 1000


# ============================================================
# Aggregate Statistics
# ============================================================

subsystem_trials = 0

equation_statistics = [[0, 0] for _ in range(mu)]

subsystem_statistics = [[0, 0] for _ in range(r)]

extra_equations = 0


# ============================================================
# Main Experimental Loop
# ============================================================

for experiment in range(1,num_experiments + 1):

    # ========================================================
    # Generate Random Ternary sLWE Instance
    # ========================================================

    A = Matrix(Zmod(M),m,n, lambda i, j: 0)

    support_sets = []

    for row in range(m):

        support = random.sample(range(n),k)

        support_sets.append(support)

        for column in support:

            A[row,column] = randint(1,M - 1)

    secret = vector(Zmod(M),[random.choice([-1, 0, 1]) for _ in range(n)])

    error = vector(Zmod(M),[random.choice([-1, 0, 1]) for _ in range(m)])

    target = (A * secret - error) 

    # ========================================================
    # Extract r Subsystems
    # ========================================================

    extracted_subsystems = []

    attempts = 0


    for _ in range(r):
        
        flag = 0
        
        while flag == 0:

            attempts += 1

            subsystem_rows = random.sample(range(m), mu)

            subsystem_equations = [A[i] for i in subsystem_rows]

            active_variables = []

            for equation in subsystem_equations:

                for position in range(n):

                    if (equation[position] != 0 and position not in active_variables):

                        active_variables.append(position)

            if len(active_variables) <= v:

                active_variables.sort()

                extracted_subsystems.append((subsystem_equations, subsystem_rows))

                flag = 1

    subsystem_trials += attempts

    # ========================================================
    # Rearrange Subsystems for Maximal Overlap
    # ========================================================

    ordering = [0]

    covered_variables = set()

    _, rows = extracted_subsystems[0]

    for row in rows:

        covered_variables.update(support_sets[row])

    while len(ordering) < r:

        best = -1
        best_overlap = -1

        for candidate in (set(range(r)) - set(ordering)):

            candidate_variables = set()

            _, rows = extracted_subsystems[candidate]

            for row in rows:

                candidate_variables.update(support_sets[row])

            overlap = len(covered_variables & candidate_variables)

            if overlap > best_overlap:

                best_overlap = overlap

                best = candidate

        ordering.append(best)

        _, rows = (extracted_subsystems[best])

        for row in rows:

            covered_variables.update(support_sets[row])

    ordered_subsystems = [extracted_subsystems[i] for i in ordering ]

    merged_table = []
    merged_variables = []

    # ========================================================
    # Solve and Merge Subsystems
    # ========================================================

    for subsystem_index in range(r):

        # ========================================================
        # Rearrange  Equations in Subsystems for Maximal Overlap
        # ========================================================

        equations, rows = (ordered_subsystems[subsystem_index])

        # Greedy reordering: start from rows[0], then pick the
        # remaining equation with the most variable overlap each time

        equation_ordering = [rows[0]]

        covered_equation_variables = set(support_sets[rows[0]])

        remaining_rows = list(rows[1:])

        while len(equation_ordering) < mu:

            best_equation = -1

            best_equation_overlap = -1

            for candidate_row in remaining_rows:

                overlap = len(set(support_sets[candidate_row]) & covered_equation_variables)

                if overlap > best_equation_overlap:

                    best_equation_overlap = overlap

                    best_equation = candidate_row

            equation_ordering.append(best_equation)

            covered_equation_variables.update(support_sets[best_equation])

            remaining_rows.remove(best_equation)

        reordered_indices = [rows.index(equation_ordering[i]) for i in range(mu)]

        local_table = []
        local_variables = []

        for equation_index in range(mu):

            eq_pos = reordered_indices[equation_index]

            solutions, variables = solve(equations[eq_pos], target[rows[eq_pos]])

            local_table, local_variables = Merge([local_table,local_variables],[solutions,variables])

            equation_statistics[equation_index][0] += len(local_variables)

            equation_statistics[equation_index][1] += len(local_table)

        merged_table, merged_variables = Merge([local_table,local_variables], [merged_table,merged_variables])


        subsystem_statistics[subsystem_index][0] += len(merged_variables)

        subsystem_statistics[subsystem_index][1] += len(merged_table)

    # ========================================================
    # Additional Equation Recovery
    # ========================================================

    additional = 0

    while len(merged_variables) < n:

        additional += 1

        best = 0
        best_gain = 0

        for row in range(m):

            gain = len(set(support_sets[row]) - set(merged_variables))

            if gain < 10:

                if gain > best_gain:

                    best = row
                    best_gain = gain

        solutions, variables = solve(A[best],target[best])

        merged_table, merged_variables = Merge([solutions,variables], [merged_table,merged_variables])

    extra_equations += additional

    # ========================================================
    # Reporting
    # ========================================================

    print("=" * 60)

    print("Experiment",experiment,"| Success:",secret in merged_table,'\n')

    

    print("Average subsystem trials:",round(subsystem_trials/(experiment* r)),'\n')

    
    # -------------------------------------------------------- 
    # Equation-level statistics 
    # -------------------------------------------------------- 
    for i in range(mu): 
        print( "Variables after", i + 1, "equations =", round( equation_statistics[i][0] / ( experiment * r ), 2 ), "| Solutions =", round( equation_statistics[i][1] / ( experiment * r ), 2 ) ) 
        print() 
    # -------------------------------------------------------- 
    # Subsystem-level statistics 
    # -------------------------------------------------------- 
    for i in range(r): 
        print( "Variables after", i + 1, "subsystems =", round( subsystem_statistics[i][0] / experiment, 2 ), "| Solutions =", round( subsystem_statistics[i][1] / experiment, 2 ) ) 
        print() 
    print( "Average extra equations =", round( extra_equations / experiment, 2 ) )






