# ----------------------------------------
# CELL SUPPRESSION PROBLEM 
# USING BENDERS DECOMPOSITION
# (using dual formulation of subproblem)
# ----------------------------------------

## GENERAL PARAMETERS ##

# Cell data
param ncells;
set CELLS := 1..ncells;
param a {CELLS};      # Original value
param lb {CELLS};     # Lower bound
param ub {CELLS};     # Upper bound
param c {CELLS};      # Suppression cost
param is_p {CELLS};   # 1 if primary, 0 otherwise

# Primary cells data
param npcells;
set PCELLS := 1..npcells;
param p {PCELLS};     # Index of the primary cell
param plpl {PCELLS};  # Lower protection level
param pupl {PCELLS};  # Upper protection level

# Table constraints (Sparse Matrix A)
param nconstraints;
set CONSTRS := 1..nconstraints;
param b {CONSTRS}; # b = (0,0,...0)
param nnz;
set MAT_ENTRIES := 1..nnz;
param begconst {1..nconstraints+1};
param coef {MAT_ENTRIES};
param xcoef {MAT_ENTRIES};

### MASTER PROBLEM: given current restrictions, minimze total suppression cost

# Benders cuts parameters
param nCUT >= 0 integer;
param gamma {1..nCUT};
param mu {CELLS, 1..nCUT};
param delta {CELLS, 1..nCUT};
param protection_lvl_k{1..nCUT};

# Master variables: 1 if suppressed
var y {CELLS} binary;

# Objective: minimize suppression cost
minimize TotalCost: sum {i in CELLS} c[i] * y[i];

# Primary cells MUST be suppressed
subj to SuppressPrimary {k in PCELLS}: y[p[k]] = 1;

subj to BendersFeasibilityCut {k in 1..nCUT}:
	(
	sum {i in CELLS} mu[i,k]*(lb[i]-a[i])*y[i]
	+ sum {i in CELLS} delta[i,k]*(ub[i]-a[i])*y[i]
	+ gamma[k] * protection_lvl_k[k]
	)
	<=0;


### SUBPROBLEM: given current suppressed cells, are the primary cells protected?

param y_fixed {CELLS};
param p_index;        # The cell index of the primary being checked
param L_or_U;         # +1 for Upper protection, -1 for Lower protection
param protection_lvl; # The required plpl or pupl

# Dual Variables
var Lambda {CONSTRS};			# Dual of Ax=0
var Mu {CELLS} >= 0;			# Dual of x >= (lb-a)y
var Delta {CELLS} <= 0;			# Dual of x <= (ub-a)y
var Gamma >= 0;					# Dual of the protection level constraint

# Objective function
maximize Dual_Objective_Function:
	sum{i in CONSTRS} Lambda[i]* b[i]
    + sum {i in CELLS} Mu[i]*(lb[i]-a[i])*y_fixed[i]
    + sum {i in CELLS} Delta[i]*(ub[i]-a[i])*y_fixed[i]
    + Gamma * protection_lvl;

# Dual Constraints (Sum of duals for each cell i)
subj to Dual_Restriction {i in CELLS}:
    (sum {l in MAT_ENTRIES: xcoef[l] == i} coef[l] * (
        sum {r in CONSTRS: begconst[r] <= l < begconst[r+1]} Lambda[r]
    ))
    + Mu[i] + Delta[i] 
    + (if i == p_index then Gamma * L_or_U else 0) == 0;