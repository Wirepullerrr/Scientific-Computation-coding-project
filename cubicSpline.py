# Natural cubic spline interpolation
import numpy as np
from numpy import size
import matplotlib.pyplot as pyp

# Returns the maximum size of x
def input_sz(x):
    _dim = np.array(x).ndim
    if _dim > 1:
        # n = max([NUM_ROWS, NUM_COLS])
        return max(size(x,0),size(x,1))
    elif _dim == 1:
        return size(x)
    else:
        return 0

# Plot cubic spline as well as its derivative and curvature
def cubicSpline(x, y,plot=True):
    n = input_sz(x)

    coefs = cubicSplineSetup(x, y)
    
    m = 20*n
    minx = min(x)
    maxx = max(x)
    step = (maxx-minx)/m
    t = np.arange(minx,maxx+step,step)

    val = np.zeros(m+1)
    der = np.zeros(m+1)
    cur = np.zeros(m+1)
    for j in range(m+1):
        [val[j], der[j], cur[j]] = cubicSplineEval(t[j], x, y, coefs)
    # just including a basic pause (there are other methods)
    if plot:
        input("Press Enter to continue . . .")
        pyp.figure()
        pyp.plot(t,val,label="spline")
        pyp.plot(x,y,'k*')
        # pyp.savefig('CubicSplinePlot.png') # save just spline
        pyp.plot(t[1:m-1],der[1:m-1],label="derivative")
        pyp.plot(t,cur,label="curvature")
        # pyp.savefig('CubicSplineDerivatives.png') # save all curves
        pyp.xlim(minx,maxx)
        pyp.legend()
        pyp.show()
    
    
# Evaluate cubic spline interpolant
def cubicSplineEval(t, x, y, coefs):
    n = input_sz(x)

    # Take care of extrapolation
    if (t <= x[0]):
        # value, deriv, curv
        return [y[0],0,0]
    if (t >= x[-1]):
        # value, deriv, curv
        return [y[-1],0,0]

    # Find interpolation interval
    i = 0 
    while (t >= x[i+1]):
        i = i+1
    
    tt = t - x[i]
    value = y[i] + tt*(coefs[0,i] + tt*(coefs[1,i]+tt*coefs[2,i]))
    deriv = coefs[0,i] + tt*(2*coefs[1,i] + tt*3*coefs[2,i])
    curv = 2*coefs[1,i] + tt*6*coefs[2,i]
    return [value,deriv,curv]


# Set up cubic spline coefficients
def cubicSplineSetup(x, y):
    n = input_sz(x)
    coefs = np.zeros(shape=(3,n))
    
    # set up tridiagonal matrix and rhs
    superDiag = np.zeros(n)
    diag = np.zeros(n)
    subDiag = np.zeros(n)
    rhs = np.zeros(n)
    diag[0] = 1
    diag[-1] = 1
    for i in range(1,n-1):
        superDiag[i] = x[i+1] - x[i]
        subDiag[i] = x[i] - x[i-1]
        diag[i] = 2*(subDiag[i] + superDiag[i])
        rhs[i] = 3*((y[i+1]-y[i])/superDiag[i] -(y[i] -y[i-1])/subDiag[i])
    print("superDiag =",superDiag)
    print("diag =",diag)
    print("subDiag =",subDiag)
    print("rhs =",rhs)

    # Solve tridiagonal system for coefs(1,i)
    for i in range(1,n):
        multiplier = -1*subDiag[i]/diag[i-1]
        diag[i] = diag[i] + multiplier*superDiag[i-1]
        rhs[i] = rhs[i] + multiplier*rhs[i-1]
    # set last coeff in second row to be 0
    coefs[1,-1] = 0
    for i in range(n-2,0,-1):
        coefs[1,i] = (rhs[i] - superDiag[i]*coefs[1,i+1])/diag[i]
    coefs[1,0] = 0

    # solve for remaining coefficients
    for i in range(n-1):
        delta = x[i+1] - x[i]
        coefs[0,i] = (y[i+1]-y[i])/delta - delta*(2*coefs[1,i] + coefs[1,i+1])/3
        coefs[2,i] = (coefs[1,i+1]-coefs[1,i])/(3*delta)
    print("Spline Coefficients")
    print("b =",coefs[0,:-1])
    print("c =",coefs[1,:-1])
    print("d =",coefs[2,:-1])
    return coefs

#### THE FOLLOWING SHOWS BASIC USAGE
##  x = [0,1,2]
##  y = [3,-2,1]
##  cubicSpline(x,y)
#### to turn off plotting 
##  cubicSpline(x,y,False)