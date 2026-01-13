### SETS AND PARAMETERS ###
param ncells;					#number of cells
param npcells;					#number of primary/sensitive cells
param nconstraints;				#number of table constraints, number of rows of matrix A
param nnz;						#number of nonzeros in matrix A

set CELLS := 1..ncells;			# Index set of all cells
set P_INDEX := 1..npcells; 		# Index set of primary/sensitive cells
set CONSTR := 1..nconstraints;	# Index set of table consistency constraints
set NZ := 1..nnz;				# Index set of nonzero coefficients in sparse matrix A

param a {CELLS};				#original values of cells
param lb {CELLS};				#lower bounds of cells' value
param ub {CELLS};				#upper bounds of cells' value
param c {CELLS};				#cost of suppressing cells
param is_p {CELLS};				#whether cell is primary/sensitive cell

param p {P_INDEX};      		# Maps primary index to cell number
param plpl {P_INDEX};  		 	# Lower protection level of primary cells
param pupl {P_INDEX};   		# Upper protection level of primary cells

#To build Ax= b table constraints
param b {CONSTR};				# right-hand side vector of table constraints 
param begconst {1..nconstraints+1}; #pointer array for sparse A
param coef {NZ};				#values of the nonzero coefficients in A
param xcoef {NZ};				#column indices of the nonzero coefficients in A



# --- Benders Cut Storage ---
param nCUT >= 0 integer;		#number of cuts created by subproblem

# dual multipliers calculated by subproblem
param mu_l {CELLS, 1..nCUT};	# dual multipliers for lower protection of cells for each cut 
param mu_u {CELLS, 1..nCUT};	# dual multiplier for upper constraints of cells for each cut 
param protection_lvl {1..nCUT};	# protection level of each cut 
param gamma {1..nCUT};			# dual multiplier for protection level for each cut


var y {CELLS} binary;                # cells to suppress
# ---------------------------------------------------------
# MASTER PROBLEM
# ---------------------------------------------------------

# minimize total suppression cost
minimize Total_Cost: sum {i in CELLS} c[i] * y[i];

# each primary cell must be suppressed
subject to Primary_Must_Be_Suppressed {k in P_INDEX}:
    y[p[k]] = 1;
    
# ensures protection feasibility for all previously detected violations
# the benders feasibility cuts are generated from the subproblem dual unbounded rays
subj to ProtectionCut {k in 1..nCUT}:
    sum {i in CELLS} (mu_l[i,k]*(lb[i]-a[i]) 
    + mu_u[i,k]*(ub[i]-a[i])) * y[i]
    + gamma[k] * protection_lvl[k]
    <= 0; 

# ---------------------------------------------------------
# DUAL SUBPROBLEM 
# ---------------------------------------------------------

param y_fixed {CELLS}; # cells suppressed calculated by master 
param p_index;         # cell index of the primary cell being checked
param L_or_U;          # decides sign of protection level depending on cut_type
						#-1 for Lower protection (min xp), +1 for Upper (max xp)
param PROTECTION_LVL;	# protection level of primary cell being checked

var LAMBDA {CONSTR};               # Dual variables for right-hand side of table constraints
var MU_L {CELLS} >= 0;             # Dual variables for lower protection
var MU_U {CELLS} <= 0;             # Dual variables for upper protection
var GAMMA >= 0;						# Dual variable for protection level

# unboundedness indicates violation of protection requirements
maximize Dual_Objective_Function:
    sum {i in CONSTR} LAMBDA[i]* b[i]
    + sum {i in CELLS} MU_L[i]*(lb[i]-a[i])*y_fixed[i]
    + sum {i in CELLS} MU_U[i]*(ub[i]-a[i])*y_fixed[i]
    + GAMMA * PROTECTION_LVL;

# dual feasibility constraints linking table consistency,
# cell bounds, and protection requirement for the current primary cell
subj to Dual_Restriction {i in CELLS}:
    (sum {l in NZ: xcoef[l] == i} coef[l] * (
        sum {r in CONSTR: begconst[r] <= l < begconst[r+1]} LAMBDA[r]
    ))
    + MU_L[i] + MU_U[i] 
    + (if i == p_index then GAMMA * L_or_U else 0) == 0;
   