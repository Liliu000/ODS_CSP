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
example: csp problem data examples
csp_benders.mod: contains subproblem and master problem (needs to specify bender cuts-doesnot run)
csp_benders.run: script to implement solution
ONLY DEV: 
-   csp_sol_gemini.txt: solution to example_2D of gemini
-   Advice teacher: First tries to debug ample use unitary cost 

# Git warnings??
warning: in the working copy of 'examples/example_2D.ampl', LF will be replaced by CRLF the next time Git touches it
warning: in the working copy of 'examples/targus.ampl', LF will be replaced by CRLF the next time Git touches it