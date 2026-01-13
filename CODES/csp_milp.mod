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

### VARIABLES ###
var y {CELLS} binary;
var x_l {CELLS, P_INDEX};
var x_u {CELLS, P_INDEX};

param y_fixed {CELLS}; #for aftwer when auditing

### OBJECTIVE ###
minimize Total_Cost:
    sum {i in CELLS} c[i] * y[i];

### CONSTRAINTS ###
subject to Lower_Bound_Prot {j in P_INDEX}:
    x_l[p[j], j] <= a[p[j]] - plpl[j]; # Modified to match image interval

subject to Upper_Bound_Prot {j in P_INDEX}:
    x_u[p[j], j] >= a[p[j]] + pupl[j]; # Modified to match image interval

subject to Flow_Conservation_L {k in CONSTR, j in P_INDEX}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x_l[xcoef[l], j] = b[k];

subject to Flow_Conservation_U {k in CONSTR, j in P_INDEX}:
    sum {l in begconst[k]..begconst[k+1]-1} coef[l] * x_u[xcoef[l], j] = b[k];

subject to Capacity_L_Lower {i in CELLS, j in P_INDEX}:
    a[i] + (lb[i] - a[i]) * y[i] <= x_l[i, j];

subject to Capacity_L_Upper {i in CELLS, j in P_INDEX}:
    x_l[i, j] <= a[i] + (ub[i] - a[i]) * y[i];

subject to Capacity_U_Lower {i in CELLS, j in P_INDEX}:
    a[i] + (lb[i] - a[i]) * y[i] <= x_u[i, j];

subject to Capacity_U_Upper {i in CELLS, j in P_INDEX}:
    x_u[i, j] <= a[i] + (ub[i] - a[i]) * y[i];

subject to Primary_Must_Be_Suppressed {i in CELLS: is_p[i] == 1}:
    y[i] = 1;