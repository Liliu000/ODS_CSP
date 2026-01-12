### SETS AND PARAMETERS ###
param ncells;
param npcells;
param nconstraints;
param nnz;

set CELLS := 1..ncells;
set P_INDEX := 1..npcells; 
set CONSTR := 1..nconstraints;
set NZ := 1..nnz;

param a {CELLS};
param lb {CELLS};
param ub {CELLS};
param c {CELLS};
param is_p {CELLS};

param p {P_INDEX};      # Maps primary index to cell number
param plpl {P_INDEX};   # Lower protection level
param pupl {P_INDEX};   # Upper protection level

param b {CONSTR};
param begconst {1..nconstraints+1};
param coef {NZ};
param xcoef {NZ};



# --- Benders Cut Storage ---
param nCUT >= 0 integer;
param mu_l {CELLS, 1..nCUT};
param mu_u {CELLS, 1..nCUT};
param protection_lvl {1..nCUT};
param gamma {1..nCUT};


var y {CELLS} binary;                # Master variables
# ---------------------------------------------------------
# MASTER PROBLEM
# ---------------------------------------------------------
minimize Total_Cost: sum {i in CELLS} c[i] * y[i];

subject to Primary_Must_Be_Suppressed {k in P_INDEX}:
    y[p[k]] = 1;
    
subj to ProtectionCut {k in 1..nCUT}:
    sum {i in CELLS} (mu_l[i,k]*(lb[i]-a[i]) 
    + mu_u[i,k]*(ub[i]-a[i])) * y[i]
    + gamma[k] * protection_lvl[k]
    <= 0; 


# ---------------------------------------------------------
# DUAL SUBPROBLEM 
# ---------------------------------------------------------

param y_fixed {CELLS};
param p_index;        # The cell index of the primary being checked
param L_or_U;         # -1 for Lower protection (min xp), +1 for Upper (max xp)
param PROTECTION_LVL;

var LAMBDA {CONSTR};               # Dual variables for lower protection
var MU_L {CELLS} >= 0;             # Dual variables for lower protection
var MU_U {CELLS} <= 0;             # Dual variables for lower protection
var GAMMA >= 0;

maximize Dual_Objective_Function:
    sum {i in CONSTR} LAMBDA[i]* b[i]
    + sum {i in CELLS} MU_L[i]*(lb[i]-a[i])*y_fixed[i]
    + sum {i in CELLS} MU_U[i]*(ub[i]-a[i])*y_fixed[i]
    + GAMMA * PROTECTION_LVL;

# Dual Restriction: Ensures dual feasibility relative to cell deviations [cite: 584, 596]
subj to Dual_Restriction {i in CELLS}:
    (sum {l in NZ: xcoef[l] == i} coef[l] * (
        sum {r in CONSTR: begconst[r] <= l < begconst[r+1]} LAMBDA[r]
    ))
    + MU_L[i] + MU_U[i] 
    + (if i == p_index then GAMMA * L_or_U else 0) == 0;
   