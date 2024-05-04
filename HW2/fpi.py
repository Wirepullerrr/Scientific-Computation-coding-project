#!/usr/bin/env python3
'''
 FIXED POINT (PICARD) ITERATION METHOD

 Solves the problem
   g(x) = x
 using fixed point iteration. For a known true solution calculates errors

 The main function is fpi:
 
  [state, x, errors, iter] = fpi(@g, x0, tolerance, maxIteration, debug);

  Inputs:
    @g            Handle to function g
    x0            The initial guess at the fixed point
    tolerance     The convergence tolerance (must be > 0).
    maxIteration  The maximum number of iterations that can be taken.
    debug         Boolean for printing out information on every iteration.
  Outputs:
    x             The solution
    errors        Array with errors at each iteration
    iter          number of iterations to convergence
  Return:
    state         An error status code.
      SUCCESS     Sucessful termination.
      WONT_STOP   Error: Exceeded maximum number of iterations.
'''
from math import sqrt
from numpy import zeros
from enum import Enum
# NOTE: You must have the Enum package. If you don't , please use fpi2.py
############################## VARIABLES #############################
class STATE(Enum):
    SUCCESS = 0
    WONT_STOP = 1

x_true = sqrt(2)
############################## FUNCTIONS #############################
#f = lambda x:(x-0.5*(x^2-2))
def f(x):
    return x - 0.5*(x*x-2)

def fpi(func,x0,TOL,MAX_ITERS,debug):
    global x_true
    prec = 12
    fmt = f"Iter %d: x= %.{prec}g, error = %.{prec}g"

    errors = zeros(MAX_ITERS+1)
    x = x0
    err = abs(x - x_true)
    errors[0] = err
    
    if debug:
        print(fmt%(0,x,err))
    
    ## FPI Loop
    for itn in range(1,MAX_ITERS+1):
        gx = func(x)
        dx = abs(x-gx)
        x = gx
        err = abs(x - x_true)
        errors[itn] = err
        if debug:
            print(fmt % (itn, x, err))
        
        # Check error tolerance
        if (dx <= TOL*(abs(x)+1)):
            iter = itn
            state = STATE.SUCCESS
            return state,x,errors, iter
    
    state = STATE.WONT_STOP
    iter = itn
    return state,x, errors, iter


################################ MAIN ###############################

###input
print("Solve the problem g(x)=x using fixed point iteration")
x0 = float(input("Enter guess at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve for the fixed point
[s,x,errors,iters] = fpi(f,x0,tol,maxIter,debug)
if s is STATE.SUCCESS:
    print("The root is %.16g"%(x))
    print("f(%.16g) = %.16g"%(x,f(x)))
    print("The number of iterations is %d"%(iters))
elif s is STATE.WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
else:
    print("ERROR: Coding error!")
    exit(1)

errors = errors[:iters+1] ## double-check this
input("Press ENTER to continue . . .") # basically just adds a pause
print("errors =",errors)
