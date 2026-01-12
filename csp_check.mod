# Parameters to be filled by the run script
param target_p;
#param y_fixed {CELLS};

# The variable the attacker tries to manipulate
var x_audit {CELLS};

# Objectives from Image 2
minimize Min_Cell_Val: x_audit[target_p];
maximize Max_Cell_Val: x_audit[target_p];

# Constraints based on Image 2
subject to Audit_Conservation {k in CONSTR}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x_audit[xcoef[l]] = b[k];

subject to Suppressed_Bounds {i in CELLS: y_fixed[i] > 0.5}:
    lb[i] <= x_audit[i] <= ub[i];

subject to Unsuppressed_Fixed {i in CELLS: y_fixed[i] <= 0.5}:
    x_audit[i] = a[i];