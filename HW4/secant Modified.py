#!/usr/bin/env python3
"""
 SECANT METHOD

 Solves the problem
    f(x) = 0
 using Secant method. For a known true solution calculates errors.

 The main function is secant:

 [state, x, errors] = secant(g0, g1,tolerance, maxIteration, debug)

 Inputs:
   g0             The initial guess at the solution.
   g1              The second guess at the solution.
   tolerance      The convergence tolerance (must be > 0).
   maxIteration   The maximum number of iterations that can be taken.
   debug          Boolean to set debugging output.
 Outputs:
   x              The solution.
 Return:
   state          An error status code.
     SUCCESS      Successful termination.
     WONT_STOP    Error: Exceeded maximum number of iterations.
     BAD_ITERATE  Error: The function had a vanishing derivative.

  Remark: We assume we are given the function
   f              The name of the function for which a root is sought.
"""
from math import sqrt
from numpy import zeros
############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_ITERATE = 2
x_true = 0.20518292
############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4
        #4000 - 20000 * (x * (1 + x) ** 6 / ((1 + x) ** 6 - 1)) AD4

def secant(g0, g1,TOL,MAX_ITERS,debug):
    global x_true, SUCCESS, WONT_STOP, BAD_ITERATE
    prec = 12
    eps = 1e-20
    # formatting string, this decides how output will look
    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g"

    errors = zeros(MAX_ITERS+2)
    x = g1
    f0 = f(g0)
    errors[0] = abs(g0-x_true)
    errors[1] = abs(g1-x_true)
    ratios_l = zeros(MAX_ITERS)
    ratios_sl = zeros(MAX_ITERS)
    
    if debug:
        print("Guess 0: x=%f, error=%.8g"%(g0,errors[0]))
        print("Guess 1: x=%f, error=%.8g"%(g1,errors[1]))
    
    ## Secant Loop
    for itn in range(1,MAX_ITERS+1):
        fx = f(x)
        if(abs(fx-f0) < eps):
            state = BAD_ITERATE
            iters = itn
            return state,x,errors,iters,ratios_l[:itn-1], ratios_sl[:itn-2]

        dx = -f(x)*(x-g0)/(fx-f0)
        g0 = x
        f0 = fx
        x += dx
        err = abs(x - x_true)
        errors[itn+1] = err

        if itn > 1:
            ratios_l[itn - 1] = errors[itn + 1]/ errors[itn]
            if itn > 2:
                ratios_sl[itn - 2] = errors[itn+1] / (errors[itn] * errors[itn - 1])
        if debug:
            print(fmt % (itn, x, dx,err))
        
        # Check error tolerance
        if (abs(dx) <= TOL):
            iter = itn
            state = SUCCESS
            return state,x,errors, iter, ratios_l[:itn], ratios_sl[:itn-1]
    
    state = WONT_STOP
    iter = itn
    return state,x, errors, MAX_ITERS, ratios_l[:MAX_ITERS], ratios_sl[:MAX_ITERS-1]


################################ MAIN ###############################
###input
print("Solve the problem f(x)=0 using Secant method")
x0 = float(input("Enter guess 0 at root: "))
x1 = float(input("Enter guess 1 at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve 
[s,x,errors,iters, ratios_l, ratios_sl] = secant(x0,x1,tol,maxIter,debug)
if s == SUCCESS:
    print(f"The root is {x:.6f}")
    print("The number of iterations is %d"%(iters))
elif s == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s == BAD_ITERATE:
    print("ERROR: Obtained a vanishing derivative!")
    exit(1)
else:
    print("ERROR: Coding error!")
    exit(1)
    
errors = errors[:iters+2]
print("errors =",errors)
print("Linear error ratios (r_l):", ratios_l)
print("Superlinear error ratios (r_sl):", ratios_sl)
