#!/usr/bin/env python3
'''
FIXED POINT (PICARD) ITERATION METHOD

Solves the problem g(x) = x using fixed point iteration.

The main function is fpi:

[state, x, errors, iter, error_ratios] = fpi(g, x0, tolerance, maxIteration, debug);

Inputs:
  g              Handle to function g
  x0             The initial guess at the fixed point
  tolerance      The convergence tolerance (must be > 0).
  maxIteration   The maximum number of iterations that can be taken.
  debug          Boolean for printing out information on every iteration.
Outputs:
  x              The solution
  errors         Array with errors at each iteration
  iter           number of iterations to convergence
  error_ratios   Array with error ratios e_i/e_{i-1} for each iteration
Return:
  state          An error status code.
    SUCCESS      Successful termination.
    WONT_STOP    Error: Exceeded maximum number of iterations.
'''
import numpy as np

SUCCESS = 0
WONT_STOP = 1


# Define g(x) based on the given equation 1 - 6x^3 = e^{2x} - 5x
# This needs to be adjusted to an appropriate form for fixed point iteration
def g(x):
    # Placeholder function, adjust according to your derivation for the equation
    return np.cbrt((-np.exp(2 * x) + 5 * x + 1) / 6)
    #(np.exp(2*x)-1+6*x**3)/5 this is function2

def fpi(func, x0, TOL, MAX_ITERS, debug):
    errors = np.zeros(MAX_ITERS + 1)
    error_ratios = np.zeros(MAX_ITERS)  # To store error ratios
    x = x0
    errors[0] = np.abs(func(x0) - x0)  # Initial error

    if debug:
        print(f"Iter 0: x= {x:.6f}, error = {errors[0]:.6f}")

    for itn in range(1, MAX_ITERS + 1):
        gx = func(x)
        err = np.abs(gx - x)
        errors[itn] = err
        error_ratios[itn - 1] = errors[itn] / errors[itn - 1] if itn > 1 else np.inf

        if debug:
            print(f"Iter {itn}: x= {gx:.6f}, error = {err:.6f}, error ratio = {error_ratios[itn - 1]:.6f}")

        if err <= TOL:
            return SUCCESS, gx, errors[:itn + 1], itn, error_ratios[:itn]

        x = gx

    return WONT_STOP, x, errors[:itn + 1], itn, error_ratios[:itn]


### MAIN
print("Solve the problem g(x)=x using fixed point iteration")
x0 = float(input("Enter guess at root: "))
tol = float(input("Enter tolerance: "))
maxIter = int(input("Enter maxIteration: "))
debug = input("Monitor iterations? (1/0): ") == '1'

state, x, errors, iters, error_ratios = fpi(g, x0, tol, maxIter, debug)
if state == SUCCESS:
    print(f"The root is {x:.16g}")
    print(f"The number of iterations is {iters}")
else:
    print(f"ERROR: Failed to converge in {maxIter} iterations!")

# Optionally, print errors and error ratios if needed
print("errors =", errors)
print("error ratios =", error_ratios)




#additional problem 3
from numpy import log
import matplotlib.pyplot as pyp
x = log(errors[:-1])
y = log(errors[1:])
dx = x[1:]-x[:-1]
dy = y[1:]-y[:-1]
slopes = [dy[i]/dx[i] for i in range(len(dx))]
print("slopes = ",slopes)
pyp.plot(x,y,"bo-")
pyp.xlabel("log(e_i)")
pyp.ylabel("log(e_{i+1})")
pyp.grid(True)
# This saves to a file
pyp.savefig("./LogErrorsPlot.png")
# This shows it on your screen
pyp.show()


