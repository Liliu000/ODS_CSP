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

# Benders Cut 
set L_CUT_IDS default {}; 
param l_mu_l {L_CUT_IDS, CELLS};
param l_mu_u {L_CUT_IDS, CELLS};
param l_cut_p_idx {L_CUT_IDS};

set U_CUT_IDS default {};
param u_mu_l {U_CUT_IDS, CELLS};
param u_mu_u {U_CUT_IDS, CELLS};
param u_cut_p_idx {U_CUT_IDS};

#DUAL
param y_fixed {CELLS};
param target_p_cell_idx; # The cell number p[j] being checked

# VARIABLES
# Master variables
var y {CELLS} binary;                
# Dual variables for lower
var lambda_L {CONSTR};              
var mu_l_L {CELLS} >= 0;             
var mu_u_L {CELLS} >= 0;             
# Dual variables for upper
var lambda_U {CONSTR};               
var mu_l_U {CELLS} >= 0;             
var mu_u_U {CELLS} >= 0;             

# MASTER PROBLEM
minimize Total_Cost: sum {i in CELLS} c[i] * y[i];

subject to Primary_Must_Be_Suppressed {i in CELLS: is_p[i] == 1}:
    y[i] = 1;

#CUTS
subject to Lower_Benders_Cuts {k in L_CUT_IDS}:
    -plpl[l_cut_p_idx[k]] >= sum {i in CELLS} (
        (lb[i] - a[i]) * l_mu_l[k, i] - (ub[i] - a[i]) * l_mu_u[k, i]
    ) * y[i];

subject to Upper_Benders_Cuts {k in U_CUT_IDS}:
    pupl[u_cut_p_idx[k]] <= sum {i in CELLS} (
        -(lb[i] - a[i]) * u_mu_l[k, i] + (ub[i] - a[i]) * u_mu_u[k, i]
    ) * y[i];

# DUAL SUBPROBLEM LOWER 
maximize Dual_Obj_L: 
    sum {i in CELLS} ((lb[i] - a[i]) * mu_l_L[i] - (ub[i] - a[i]) * mu_u_L[i]) * y_fixed[i];

# A^T * lambda + mu_l - mu_u = e_p
subject to Dual_Constraint_L {i in CELLS}:
    sum {k in CONSTR, l in begconst[k]..begconst[k+1]-1: xcoef[l] == i} (coef[l] * lambda_L[k]) 
    + mu_l_L[i] - mu_u_L[i] = (if i = target_p_cell_idx then 1 else 0);

# DUAL SUBPROBLEM UPPER  

maximize Dual_Obj_U: 
    sum {i in CELLS} (-(lb[i] - a[i]) * mu_l_U[i] + (ub[i] - a[i]) * mu_u_U[i]) * y_fixed[i];

# A^T * lambda + mu_l - mu_u = -e_p
subject to Dual_Constraint_U {i in CELLS}:
    sum {k in CONSTR, l in begconst[k]..begconst[k+1]-1: xcoef[l] == i} (coef[l] * lambda_U[k]) 
    + mu_l_U[i] - mu_u_U[i] = (if i = target_p_cell_idx then -1 else 0); #<--do not understand ep