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
param pupl {P_INDEX};   # Upper protection level (corrected name from image)

param b {CONSTR};
param begconst {1..nconstraints+1};
param coef {NZ};
param xcoef {NZ};

### VARIABLES ###
var y {CELLS} binary;

# Flow variables for lower (l) and upper (u) scenarios for each primary p
var x_l {CELLS, P_INDEX};
var x_u {CELLS, P_INDEX};

### OBJECTIVE ###
minimize Total_Cost:
    sum {i in CELLS} c[i] * y[i];

### CONSTRAINTS ###

# 1. Protection Level Constraints for each Primary Cell p
# x_p <= -lpl_p and x_p >= upl_p
subject to Lower_Bound_Prot {j in P_INDEX}:
    x_l[p[j], j] <= -plpl[j];

subject to Upper_Bound_Prot {j in P_INDEX}:
    x_u[p[j], j] >= pupl[j];

# 2. Conservation Constraints (Ax = 0) using sparse matrix indexing
subject to Flow_Conservation_L {k in CONSTR, j in P_INDEX}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x_l[xcoef[l], j] = b[k];

subject to Flow_Conservation_U {k in CONSTR, j in P_INDEX}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x_u[xcoef[l], j] = b[k];

# 3. Capacity Constraints linking x and y
# (l_i - a_i) * y_i <= x_i <= (u_i - a_i) * y_i
subject to Capacity_L_Lower {i in CELLS, j in P_INDEX}:
    (lb[i] - a[i]) * y[i] <= x_l[i, j];

subject to Capacity_L_Upper {i in CELLS, j in P_INDEX}:
    x_l[i, j] <= (ub[i] - a[i]) * y[i];

subject to Capacity_U_Lower {i in CELLS, j in P_INDEX}:
    (lb[i] - a[i]) * y[i] <= x_u[i, j];

subject to Capacity_U_Upper {i in CELLS, j in P_INDEX}:
    x_u[i, j] <= (ub[i] - a[i]) * y[i];

# 4. Mandatory suppression of primary cells
subject to Primary_Must_Be_Suppressed {i in CELLS: is_p[i] == 1}:
    y[i] = 1;