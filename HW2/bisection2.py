#!/usr/bin/env python3
'''
 BISECTION METHOD

 Solves the problem
   f(x) = 0
 using the bisection algorithm.

 [state,root] = bisection(a, b, tolerance, maxIteration, debug)

 Inputs:
   a,b           The initial bounding interval, with a root between.
   tolerance     The convergence tolerance (must be > 0).
   maxIteration  The maximum number of iterations that can be taken.
   debug         Boolean for printing out information on every iteration.
 Outputs:
   root             The solution.
   state         An error status code.
     SUCCESS     Sucessful termination.
     WONT_STOP   Error: Exceeded maximum number of iterations.
     BAD_DATA    Error: The interval may not bracket a root.

 We assume the function f is given:
   double f(double);
'''

from math import exp,cos
############################## VARIABLES #############################
SUCCESS = 0
WONT_STOP = 1
BAD_DATA = 2
BAD_ITERATE = 3
############################## FUNCTIONS #############################

def f(x):
    # return cos(x) - x
    # return x*x - 2
    return x - exp(-x)

def sgn(x):
    if (x>0):
        return 1
    elif (x<0):
        return -1
    elif (x==0):
        return 0
    else:
        return x #x is not a number

# Returns two values, the first is the error state and the second is the found root
def bisection(a,b,tolerance, maxIteration,debug):
    global SUCCESS, WONT_STOP, BAD_DATA
    x = None
    # format string
    prec = 8
    fmt = f"Iter %d: x= %.{prec}g, dx = %.{prec}g, a = %.{prec}g, b = %.{prec}g, f(x) = %.{prec+4}g"
    
    # if necessary, swap a and b
    if ( a > b):
        c = a
        a = b
        b = c
    
    fa = f(a)
    fb = f(b)

    # make sure there is a root between a and b
    if(sgn(fa)*sgn(fb)>0.0):
        return BAD_DATA,x
    
    # iteration loop
    dx = b-a
    for iteration in range(maxIteration):
        dx/=2
        x = a+dx
        fx = f(x)
        if debug:
            print(fmt % (iteration, x, dx,a,b,fx))
        
        # Check error tolerance
        if (dx <= tolerance):
            return SUCCESS,x
        
        if(sgn(fa)*sgn(fx)>0.0):
            a = x
            fa = fx
        else:
            b = x
            fb = fx
    return WONT_STOP,x
################################ MAIN ###############################

### input
print("Solves the problem f(x) = 0 on interval [a,b] using the bisection algorithm")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
tol = float(input("Enter tolerance: "))
maxIter = float(input("Enter maxIteration: "))
debug = int(input("Monitor iterations? (1/0): "))
if debug != 0: debug = 1 # technically not necessary but just some error handling

### Solve for a root
[s,root] = bisection(a,b,tol,maxIter,debug)
if s == SUCCESS:
    print("The root is %.12g"%(root))
    print("f(%.12g) = %.12g"%(root,f(root)))
    exit()
elif s == WONT_STOP:
    print("ERROR: Failed to converge in %d iterations!"%(maxIter))
elif s == BAD_DATA:
    print("ERROR: Unsuitable interval!")
else:
    print("ERROR: Coding error!")
exit(1)
