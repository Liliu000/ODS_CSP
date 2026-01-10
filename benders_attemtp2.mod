




## SUBPROBLEM ##

# Parameters

# Variables

var U{M}

maximize 




## MASTER ##

# Parameters (fixed)
param c {CELLS};

# Parameters (from subproblem)
param sl {CELLS, P, 1..2};
param su {CELLS, P, 1..2};
param l {CELLS};
param u {CELLS};

# Variables to optimize
var Y {CELLS} binary;
var Z;


minimize Total_Cost:
    Z;

# Points
subject to Lower_Bound_p{k in 1..nCUT, i in CELLS, p in P, h in 1..2}:
    Z >= (sum {j in CELLS} c[j] * Y[j]) + 
         (if cut_type[k] = "point" then (sl[i,p,h] - (l[i] - a[i])) else 0);

subject to Upper_Bound_p{k in 1..nCUT, i in CELLS, p in P, h in 1..2}:
    Z >= (sum {j in CELLS} c[j] * Y[j]) + 
         (if cut_type[k] = "point" then (su[i,p,h] - (a[i] - u[i])) else 0);

# Rays
subject to Lower_Bound_r{k in 1..nCUT, i in CELLS, p in P, h in 1..2}:
    if cut_type[k] = "ray" then
        sl[i,p,h] - (l[i] - a[i]) <= 0

subject to Upper_Bound_r{k in 1..nCUT, i in CELLS, p in P, h in 1..2}:
    if cut_type[k] = "ray" then
        su[i,p,h] - (a[i] - u[i]) <= 0