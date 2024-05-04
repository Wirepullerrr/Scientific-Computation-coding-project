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
    df              The name of the derivative of the function. The derivative
                    of f must be computed by hand and coded correctly as df, or
                    newton will not work!
'''
from math import sin,cos,pi
from numpy import zeros
from enum import Enum
# NOTE: You must have the Enum package. If you don't , please use newton2.py
############################## VARIABLES #############################
class STATE(Enum):
    SUCCESS = 0
    WONT_STOP = 1
    BAD_ITERATE = 2
x_true = 0.5
############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4

# The name of th derivative of the function. The derivative of f must 
# be computed by hand and coded correctly as df, or newton will not work!
def df(x):
    return 324*x**5 + 225*x**4 - 408*x**3 - 207*x**2 + 70*x + 16

def newton(x0,TOL,MAX_ITERS,debug):
    global x_true
    prec = 12
    eps = 1e-20
    # formatting string, this decides how output will look
    fmt = f"Iter %d: x= %.{prec}g, dx= %.{prec}g, error = %.{prec}g"

    errors = zeros(MAX_ITERS+1)
    x = x0
    err = abs(x - x_true)
    errors[0] = err
    
    if debug:
        print("Guess: x=%.8g, error=%.8g"%(x,err))
    
    ## Newton Loop
    for itn in range(1,MAX_ITERS+1):
        dfx = df(x)
        if(abs(dfx) < eps):
            state = STATE.BAD_ITERATE
            iters = itn
            return state,x,errors,iters

        dx = -f(x)/dfx
        # Use for Modified Netwon, multiplicity 2
        # dx = 2*dx
        x += dx
        err = abs(x - x_true)
        errors[itn] = err
        if debug:
            print(fmt % (itn, x, dx,err))
        
        # Check error tolerance
        if (abs(dx) <= TOL):
            iter = itn
            state = STATE.SUCCESS
            return state,x,errors, iter
    
    state = STATE.WONT_STOP
    iter = itn
    return state,x, errors, iter


################################ MAIN ###############################

###input
print("Solve the problem f(x)=0 using Newton's method")
x0 = float(input("Enter guess at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve 
[s,x,errors,iters] = newton(x0,tol,maxIter,debug)
if s is STATE.SUCCESS:
    print("The root is %.16g"%(x))
    print("The number of iterations is %d"%(iters))
    errors = errors[:iters+1]
    print("errors =",errors)
    exit()
elif s is STATE.WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s is STATE.BAD_ITERATE:
    print("ERROR: Obtained a vanishing derivative!")
else:
    print("ERROR: Coding error!")
exit(1) #technically, not necessary but good for general practice
