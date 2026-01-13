# Statement
Implement in AMPL a Benders decomposition for the cell suppression problem (CSP problem). You may test your code with the instances in the accompanying zip file.The data for each instance is:
- ncells: number of cells
- a, lb, ub, c, is_p: cell value, lower bound, upper bound, suppression cost, is primary (1 if primary, 0 otherwise)
- npcells: number of primary/sensitive cells.
- p, plpl, pupl: p[i]=j means primary i is cell number j; lower and upper protection level of each primary.
- ncosntraints: number of table constraints.
- b: right hand side of constraints Ax=b.
- nnz: number of nonzeros in matrix A.
- begconst, coef, xcoef: packed rowwise sparse matrix A. The information of constraint i is in positions begconst[i] and begconst[i+1]-1 of coef and xcoef. Coef[l], begconst[i] ≤ l ≤ begconst[i+1]-1  gives the coefficients in A of row i;
xcoef[l],  begconst[i] ≤ l ≤ begconst[i+1]-1  gives the columns of A where the coefficients are located.

Instance example_2D.ampl contains a 5x6 bidimensional table (including total row and column), where the cells are sorted rowwise, assuming that the total row and total column are the first row and column. The other files contain two additional small instances.
 
# Structure
- CODES: folder with AMPL codes for solving CSP implementing benders (csp_benders) or MILP (csp_milp), both audit the answers via csp_check
- examples: folder with data examples 
- report: latex files used to create report

