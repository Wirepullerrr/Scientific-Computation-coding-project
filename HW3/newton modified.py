#!/usr/bin/env python3
'''
 NEWTON'S METHOD

 Solves the problem f(x)=0 using Newton's method. For a known tru solution
 calculates errors.

 The main function is newton:
 
  [state,x,errors,iters] = newton(x0, tolerance, maxIteration, debug)

  Inputs:
    x0              The initial guess at the solution
    tolerance       The convergence tolerance (must be > 0).
    maxIteration    The maximum number of iterations that can be taken.
    debug           Boolean to set debugging output
  Outputs:
    x               The solution
    errors          Array with errors at each iteration
    iter            number of iterations to convergence
  Return:
    state           An error status code.
      SUCCESS       Sucessful termination.
      WONT_STOP     Error: Exceeded maximum number of iterations.
      BAD_ITERATE   Error: The function had a vanishing derivative

  Remark: We assume that we known the two functions
    f               The name of the function for which a root is sought
    df              The name of th derivative of the function. The derivative
                    of f must be computed by hand and coded correctly as df, or
                    newton will not work!
'''
from numpy import zeros, abs
############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_ITERATE = 2
x_true = -2/3
############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4

# The name of th derivative of the function. The derivative of f must 
# be computed by hand and coded correctly as df, or newton will not work!
def df(x):
    return 324*x**5 + 225*x**4 - 408*x**3 - 207*x**2 + 70*x + 16

def newton(x0,TOL,MAX_ITERS,debug):
    global x_true, SUCCESS, WONT_STOP, BAD_ITERATE
    prec = 12
    eps = 1e-20
    # formatting string, this decides how output will look
    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g, r_l = %.{prec}g, r_q = %.{prec}g"

    errors = zeros(MAX_ITERS+1)
    ratios_l = zeros(MAX_ITERS)
    ratios_q = zeros(MAX_ITERS)
    x = x0
    err = abs(x - x_true)
    errors[0] = err
    prev_error = err
    
    if debug:
        print("Guess: x=%.8g, error=%.8g"%(x,err))
    
    ## Newton Loop
    for itn in range(1,MAX_ITERS+1):
        dfx = df(x)
        if(abs(dfx) < eps):
            state = BAD_ITERATE
            iters = itn
            return state, x, errors[:itn], iters, ratios_l[:itn-1], ratios_q[:itn-2]

        dx = -f(x)/dfx
        # Use for Modified Netwon, multiplicity 2
        dx = 2*dx
        x += dx
        err = abs(x - x_true)
        errors[itn] = err
        
        # Check error tolerance
        if (abs(dx) <= TOL):
            iter = itn
            state = SUCCESS
            return state, x, errors[:itn+1], iter, ratios_l[:itn-1], ratios_q[:itn-2]

        if itn > 1:  # Ratios start from the second iteration
            ratios_l[itn - 1] = err / prev_error
        if itn > 2:  # Quadratic ratios require at least two previous errors
            ratios_q[itn - 2] = err / (prev_error ** 2)

        prev_error = err

        if debug:
            print(fmt % (itn, x, dx, err, ratios_l[itn - 1] if itn > 1 else 0, ratios_q[itn - 2] if itn > 2 else 0))

    state = WONT_STOP
    iter = itn
    return state, x, errors[:itn], iter, ratios_l[:itn-1], ratios_q[:itn-2]


################################ MAIN ###############################

###input
print("Solve the problem f(x)=0 using Newton's method")
x0 = float(input("Enter guess at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve 
[s,x,errors,iters, ratios_l, ratios_q] = newton(x0,tol,maxIter,debug)
if s == SUCCESS:
    print(f"The root is {x:.6f}")
    print("The number of iterations is %d"%(iters))
    errors = errors[:iters+1]
    print("errors =",errors)
    print("Linear error ratios r_l:", ratios_l)
    print("Quadratic error ratios r_q:", ratios_q)
    exit()
elif s == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s == BAD_ITERATE:
    print("ERROR: Obtained a vanishing derivative!")
else:
    print("ERROR: Coding error!")
exit(1) #technically not necessary, but good for general practice


