# ============================================
# CELL SUPPRESSION PROBLEM - BENDERS DECOMPOSITION
# Corrected model
# ============================================

# --------------------
# DATA
# --------------------
param ncells >= 1 integer;
set CELLS := 1..ncells;

param a {CELLS};
param lb {CELLS};
param ub {CELLS};
param c {CELLS};
param is_p {CELLS} binary;

param npcells >= 0 integer;
set PCELLS := 1..npcells;

param p {PCELLS};
param plpl {PCELLS};
param pupl {PCELLS};

param nconstraints >= 0 integer;
set CONSTRAINTS := 1..nconstraints;

param b {CONSTRAINTS};
param nnz >= 0 integer;
param begconst {1..nconstraints+1};
param coef {1..nnz};
param xcoef {1..nnz};

# SUBPROBLEM (DUAL)
param x_fixed {CELLS} binary;

var u_lower {PCELLS} >= 0;
var u_upper {PCELLS} >= 0;

maximize Sub_Objective:
    sum {i in PCELLS} (
        plpl[i] * u_upper[i] - pupl[i] * u_lower[i]
    );

subject to Sub_Dual {j in CELLS}:
    sum {i in PCELLS: p[i] = j} (u_upper[i] - u_lower[i]) 
    <= x_fixed[j];  # Checkout, not sure

# MASTER PROBLEM
param nCUT >= 0 integer default 0;
param cut_type {1..nCUT} symbolic within {"optimality","feasibility"};

param cut_u_lower {PCELLS,1..nCUT};
param cut_u_upper {PCELLS,1..nCUT};

var x {CELLS} binary;
var theta >= 0;

minimize Master_Objective:
    sum {j in CELLS} c[j] * x[j] + theta;

# Primary cells must be suppressed
subject to Primary_Suppressed {j in CELLS: is_p[j] = 1}:
    x[j] = 1;

# Table consistency constraints
subject to Table_Constraints {i in CONSTRAINTS}:
    sum {k in begconst[i]..begconst[i+1]-1}
        coef[k] * a[xcoef[k]] * (1 - x[xcoef[k]]) = b[i];

# Benders cuts
subject to Optimality_Cuts {k in 1..nCUT: cut_type[k] = "optimality"}:
    ??

subject to Feasibility_Cuts {k in 1..nCUT: cut_type[k] = "feasibility"}:
    ??