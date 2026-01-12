### SETS AND PARAMETERS FROM MILP ###
param ncells;
param nconstraints;
param nnz;
set CELLS := 1..ncells;
set CONSTR := 1..nconstraints;
set NZ := 1..nnz;

param a {CELLS};
param lb {CELLS};
param ub {CELLS};
param b {CONSTR};
param begconst {1..nconstraints+1};
param coef {NZ};
param xcoef {NZ};

# Input from the MILP result
param y_fixed {CELLS}; 

### AUDIT VARIABLES ###
var x {CELLS};

### AUDIT OBJECTIVES ###
minimize Min_Cell_Val: x[target_p];
maximize Max_Cell_Val: x[target_p];

### CONSTRAINTS (Based on Image 2) ###
param target_p; # The current primary cell being audited

# Ax = b
subject to Audit_Conservation {k in CONSTR}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x[xcoef[l]] = b[k];

# Bounds for suppressed cells (i in P U S)
subject to Suppressed_Bounds {i in CELLS: y_fixed[i] > 0.5}:
    lb[i] <= x[i] <= ub[i];

# Values for unsuppressed cells (i not in P U S)
subject to Unsuppressed_Fixed {i in CELLS: y_fixed[i] <= 0.5}:
    x[i] = a[i];