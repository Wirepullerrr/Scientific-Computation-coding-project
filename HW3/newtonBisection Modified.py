#!/usr/bin/env python3
'''
 NEWTON-BISECTION METHOD

 Solves the problem f(x)=0 using Newton-Bisection method. For a known true solution
 calculates errors.

 The main function is newton:
 
  [state,x,errors,iters] = newtonBIsection(x0, tolerance, maxIteration, debug)

  Inputs:
    a,b             The initial bounding interval, with a root between.
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
      BAD_DATA      Error: The interval may not bracket a root

  Remark: We assume that we known the two functions
    f               The name of the function for which a root is sought
    df              The name of the derivative of the function.
'''
import math
from numpy import zeros,sign

############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_DATA = 2
x_true = 2.403192
############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 4 * math.exp(-0.2*x) - 15* math.exp(-0.75*x)

# The name of the derivative of the function. 
def df(x):
    return 11.25 * math.exp(-0.75 *x) - 0.8 * math.exp(-0.2 *x)

def newtonBisection(a,b,TOL,MAX_ITERS,debug):
    global x_true, SUCCESS, WONT_STOP, BAD_DATA
    prec = 12
    eps = 1e-20
    # formatting string, this decides how output will look
    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g, "
    fmt += f"interval = [%.{prec}g,%.{prec}g], Newton? %d"

    # Swap a and b if necessary so a < b
    if (a>b):
        c = a
        a = b
        b = c
    fa = f(a)
    fb = f(b)

    # Make sure there is a root between a and b
    if(sign(fa)*sign(fb) > 0.0):
        state = BAD_DATA
        x = None
        errors = None
        iter = 0
        return state,x,errors,iter,[],[],[]

    errors = zeros(MAX_ITERS+1)
    x = a+(b-a)/2
    err = abs(x - x_true)
    errors[0] = err
    ratios_l = zeros(MAX_ITERS)
    ratios_q = zeros(MAX_ITERS)
    ratios_sl = zeros(MAX_ITERS)

    if debug:
        print("Interval = [%f,%f], guess x = %f, error = %f"%(a,b,x,err))
    
    fx = f(x)
    if(sign(fa)*sign(fx) > 0.0):
        a = x
    else:
        b = x

    ## NewtonBisection Loop
    for itn in range(1,MAX_ITERS+1):
        dfx = df(x)
        usedNewton = True
        if(abs(dfx) > eps):
            xNew = x - fx/dfx # Newton
            if(xNew < a or b < xNew):
                xNew = a + (b-a)/2 # Revert to Bisection
                usedNewton = False
        else:
            xNew = a + (b-a)/2;  # Revert to Bisection
            usedNewton = False

        fx = f(xNew)
        if(sign(fa)*sign(fx) > 0.0):
            a = xNew
        else:
            b = xNew

        dx = xNew - x
        x = xNew
        err = abs(x - x_true)
        errors[itn] = err
        ratios_l[itn - 1] = err / errors[itn - 1]
        ratios_q[itn - 1] = err / (errors[itn - 1] ** 2)
        if itn > 1:
            ratios_sl[itn - 1] = err/ (errors[itn - 1] * errors[itn - 2])

        if debug:
            print(fmt % (itn, x, dx,err,a,b,usedNewton))
        
        # Check error tolerance
        if (abs(dx) <= TOL):
            iter = itn
            state = SUCCESS
            return state, x, errors, iter, ratios_l[:itn], ratios_q[:itn], ratios_sl[:itn]
    
    state = WONT_STOP
    iter = itn
    return state, x, errors, MAX_ITERS, ratios_l[:MAX_ITERS], ratios_q[:MAX_ITERS], ratios_sl[:MAX_ITERS]


################################ MAIN ###############################

###input
print("Solve the problem f(x)=0 on interval [a,b] using Newton-Bisection method")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve 
s, x, errors, iters, r_l, r_q, r_sl = newtonBisection(a, b, tol, maxIter, debug)
if s == SUCCESS:
    print(f"The root is {x:.6f}.")
    print("The number of iterations is %d"%(iters))
    errors = errors[:iters+1]
    print("errors =",errors)
    if debug:  # Optionally print error ratios if debugging is enabled
        print("Linear error ratios r_l:", r_l)
        print("Quadratic error ratios r_q:", r_q)
        print("Superlinear error ratios r_sl:", r_sl)
    exit()
elif s == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s == BAD_DATA:
    print("ERROR: Unsuitable interval!")
else:
    print("ERROR: Coding error!")
 #technically, not necessary but good for general practice

