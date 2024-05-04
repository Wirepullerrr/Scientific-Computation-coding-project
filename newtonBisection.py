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
import numpy as np

############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_DATA = 2
x_true = 1
############################## FUNCTIONS #############################
# The function for which a root is sought
def f(x):
    return 54*x**6 + 45*x**5 - 102*x**4 - 69*x**3 + 35*x**2 + 16*x - 4

# The name of the derivative of the function. 
def df(x):
    return 324*x**5 + 225*x**4 - 408*x**3 - 207*x**2 + 70*x + 16

def newtonBisection(a, b, TOL, MAX_ITERS, debug):
    global x_true, SUCCESS, WONT_STOP, BAD_DATA
    errors = np.zeros(MAX_ITERS + 1)
    r_l = np.zeros(MAX_ITERS)  # Linear error ratios
    r_q = np.zeros(MAX_ITERS)  # Quadratic error ratios

    x = a + (b - a) / 2  # Initial guess
    errors[0] = abs(x - x_true)

    if debug:
        print(f"Initial guess: x={x:.6f}, error={errors[0]:.6f}")

    for itn in range(1, MAX_ITERS + 1):
        dfx = df(x)
        if abs(dfx) > 1e-20:
            dx = -f(x) / dfx
            x_new = x + dx
            if not a <= x_new <= b:
                x_new = a + (b - a) / 2  # Bisection step if Newton's step is out of bounds
                dx = x_new - x

        else:
            print("Derivative too small, reverting to bisection.")
            x_new = a + (b - a) / 2
            dx = x_new - x

        err = abs(x_new - x_true)
        errors[itn] = err

        r_l[itn - 1] = err / errors[itn - 1]
        r_q[itn - 1] = err / (errors[itn - 1] ** 2)

        x = x_new

        if debug:
            print(f"Iter {itn}: x={x:.6f}, error={err:.6e}, r_l={r_l[itn-1]:.6e}, r_q={r_q[itn-1]:.6e}")

        if err < TOL:
            return SUCCESS, x, errors[:itn + 1], itn, r_l[:itn], r_q[:itn]

    return WONT_STOP, x, errors, MAX_ITERS, r_l[:MAX_ITERS], r_q[:MAX_ITERS]

################################ MAIN ###############################

###input
print("Solve the problem f(x)=0 on interval [a,b] using Newton-Bisection method")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = bool(input("Monitor iterations? (1/0): "))

### Solve 
state, x, errors, iters, r_l, r_q = newtonBisection(a, b, tol, maxIter, debug)
if state == SUCCESS:
    print(f"The root is {x:.6f}.")
    print("The number of iterations is %d"%(iters))
    errors = errors[:iters+1]
    print("errors =",errors)
    if debug:  # Optionally print error ratios if debugging is enabled
        print("Linear error ratios r_l:", r_l)
        print("Quadratic error ratios r_q:", r_q)
    exit()
elif state == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif state == BAD_DATA:
    print("ERROR: Unsuitable interval!")
else:
    print("ERROR: Coding error!")
exit(1) #technically, not necessary but good for general practice
