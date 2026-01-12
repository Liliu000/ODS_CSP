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

param p {P_INDEX};      
param plpl {P_INDEX};   
param pupl {P_INDEX};   

param b {CONSTR};
param begconst {1..nconstraints+1};
param coef {NZ};
param xcoef {NZ};

# Pre-calculate Matrix A for dual constraints (efficiency)
set CELLS_IN_CONSTR {k in CONSTR} := 
    setof {l in begconst[k]..begconst[k+1]-1} xcoef[l];

set CONSTR_OF_CELL {i in CELLS} := 
    setof {k in CONSTR, l in begconst[k]..begconst[k+1]-1: xcoef[l] == i} k;

param A_coeff {k in CONSTR, i in CELLS_IN_CONSTR[k]} := 
    sum {l in begconst[k]..begconst[k+1]-1: xcoef[l] == i} coef[l];

### MASTER PROBLEM ###

param nCUT_L default 0;
param nCUT_U default 0;

# Storage for dual values from subproblems to build cuts
param lambda_L {CONSTR, P_INDEX, 1..1000} default 0;
param mu_l_L {CELLS, P_INDEX, 1..1000} default 0;
param mu_u_L {CELLS, P_INDEX, 1..1000} default 0;

param lambda_U {CONSTR, P_INDEX, 1..1000} default 0;
param mu_l_U {CELLS, P_INDEX, 1..1000} default 0;
param mu_u_U {CELLS, P_INDEX, 1..1000} default 0;

var y {CELLS} binary;

minimize Total_Cost:
    sum {i in CELLS} c[i] * y[i];

subject to Mandatory_Suppression {i in CELLS: is_p[i] == 1}:
    y[i] = 1;

# Benders Feasibility Cuts for Lower Protection
subject to Lower_Prot_Cut {j in P_INDEX, k in 1..nCUT_L}:
    sum {i in CELLS} ((lb[i]-a[i])*mu_l_L[i,j,k] - (ub[i]-a[i])*mu_u_L[i,j,k]) * y[i] 
    <= -plpl[j] - sum {idx in CONSTR} b[idx] * lambda_L[idx,j,k];

# Benders Feasibility Cuts for Upper Protection
subject to Upper_Prot_Cut {j in P_INDEX, k in 1..nCUT_U}:
    sum {i in CELLS} ((lb[i]-a[i])*mu_l_U[i,j,k] - (ub[i]-a[i])*mu_u_U[i,j,k]) * y[i] 
    >= pupl[j] - sum {idx in CONSTR} b[idx] * lambda_U[idx,j,k];


### SUBPROBLEM (DUAL FORMULATION) ###

param current_y {CELLS};   # Passed from Master
param current_p_cell;      # The cell index of the primary being checked

var lambda {CONSTR};
var mu_l {CELLS} >= 0;
var mu_u {CELLS} >= 0;

# Lower Scenario Dual: min x[p] => max Dual
maximize Dual_Obj_L:
    sum {k in CONSTR} b[k] * lambda[k] + 
    sum {i in CELLS} ((lb[i]-a[i]) * current_y[i] * mu_l[i] - (ub[i]-a[i]) * current_y[i] * mu_u[i]);

# Upper Scenario Dual: max x[p] => min Dual (or max -x[p])
# Here we formulate max x[p] directly as a dual
maximize Dual_Obj_U:
    sum {k in CONSTR} b[k] * lambda[k] + 
    sum {i in CELLS} ((lb[i]-a[i]) * current_y[i] * mu_l[i] - (ub[i]-a[i]) * current_y[i] * mu_u[i]);

subject to Dual_Constraint {i in CELLS}:
    sum {k in CONSTR_OF_CELL[i]} A_coeff[k,i] * lambda[k] + mu_l[i] - mu_u[i] = 
    (if i == current_p_cell then 1 else 0);