# Newton Divided Difference Interpolation

import numpy as np
import matplotlib.pyplot as pyp
import math
import argparse


############################## FUNCTIONS #############################

# Evaluate divided difference interpolant
def newtonEval(t,coefs,x):
    n = len(coefs)
    value = coefs[n-1]
    for i in range(n-2,-1,-1): # same as n-2, n-3, n-4, ..., 0
        value = value*(t-x[i]) + coefs[i]
    return value

# Set up divided difference coefficients
def newtonDDsetup(x,y):
    n = len(x)
    if (len(y) != n): 
        print("ERROR CODE 1: x and y are different sizes")
        exit(1)

    # DD level 0
    # coefs[i] = y[i] for i=0,1,2,...,n-1
    coefs = [y[i] for i in range(n)] 
    
    # DD higher levels (bottom to top, overwrite lower entries as they are finished)
    for level in range(1,n): # 1,2,3,4, ... n-1
        for i in range(n-1,level-1,-1): #n-1, n-2, ..., level
            dx = x[i] - x[i-level] 
            if (dx==0): exit(2)
            coefs[i] = (coefs[i]-coefs[i-1])/dx
    return coefs

def chebyshev_nodes(n, a, b):
    return [0.5 * (a + b) + 0.5 * (b - a) * math.cos((2*k - 1) * math.pi / (2 * n)) for k in range(1, n + 1)]

def domain_cos(x_value):
    x_revised = []
    for x in x_value:
        s = 1
        x = x % (2*math.pi)
        if x > math.pi:
            x = 2*math.pi - x
        if x > math.pi/2:
            s = -1
            x = math.pi - x
        x_revised.append((x, s))
    return x_revised

# x,y are 1d- arrays
def newtonDD(x,y,name):
    n = len(x)
    if (len(y) != n): exit(1)


    coefs = newtonDDsetup(x,y)


    m = 10*n
    minx = min(x)
    maxx = max(x)
    t = np.arange(minx,maxx+1,(maxx-minx)/m)
    val = newtonEval(t,coefs,x)

    # plot
    x1 = np.arange(-np.pi,3*np.pi, 0.1)
    y1 = [np.cos(x) for x in x1]
    # Plotting coine Graph
    pyp.plot(x1, y1, color='green')
    pyp.plot(t,val)
    for i in range(n):
        pyp.plot(x[i],y[i],'k*')
    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(-np.pi,3*np.pi)
    pyp.ylim(-1.1, 1.1)
    pyp.title("Newton DD interpolation")
    pyp.legend(["Cosx", name],loc="best")
    pyp.show()



def compute_error(x_value, coefs, nodes):
    x_revised = domain_cos(x_value)
    y_actual = [np.cos(val) for val in x_value]
    y_interpolated = [s*newtonEval(x, coefs, nodes) for (x, s) in x_revised]
    errors = [abs(actual - approx) for actual, approx in zip(y_actual, y_interpolated)]
    return errors

# Function to plot the errors
def plot_errors(x_vals, error1, error2):
    pyp.figure(figsize=(10, 6))
    pyp.plot(x_vals, error1, label='Error in Cos1(x)')
    pyp.plot(x_vals, error2, label='Error in Cos2(x)')
    pyp.xlabel('x')
    pyp.ylabel('Error')
    pyp.title('Error in Newton and Chebyshev Interpolations')
    pyp.legend()
    pyp.xlim([-np.pi, 3*np.pi])
    pyp.show()



############################## MAIN #############################
x_known = [0, np.pi / 4, np.pi / 2]
y_known = [np.cos(xi) for xi in x_known]
x_chebyshev = chebyshev_nodes(3, 0, math.pi/2)
y_chebyshev = [np.cos(xi) for xi in x_chebyshev]
'''y_known = [1, 1/math.sqrt(2),0]'''

coefs = newtonDDsetup(x_known,y_known)
coefs2 = newtonDDsetup(x_chebyshev,y_chebyshev)
x_interpolated = [-1000, -14, -4, -3, -2, -1, 1, 2, 3, 4, 14, 1000]
x_revised = domain_cos(x_interpolated)
y_actual = [np.cos(val) for val in x_interpolated]
y_interpolated = [s*newtonEval(x, coefs, x_known) for (x, s) in x_revised]
errors = [abs(actual - approx) for actual, approx in zip(y_actual, y_interpolated)]

y_interpolated2 = [s*newtonEval(x, coefs2, x_chebyshev) for (x, s) in x_revised]
errors2 = [abs(actual - approx) for actual, approx in zip(y_actual, y_interpolated2)]
print("{:>10} {:>20} {:>20} {:>20}".format("x", "cos(x)", "cos1(x)", "Error"))
for i, val in enumerate(x_interpolated):
    print("{:10.1f} {:20.6f} {:20.6f} {:20.6f}".format(val, y_actual[i], y_interpolated[i], errors[i]))

print()
print("{:>10} {:>20} {:>20} {:>20}".format("x", "cos(x)", "cos2(x)", "Error"))
for i, val in enumerate(x_interpolated):
    print("{:10.1f} {:20.6f} {:20.6f} {:20.6f}".format(val, y_actual[i], y_interpolated2[i], errors2[i]))

newtonDD(x_known, y_known, "Cos1x")
newtonDD(x_chebyshev, y_chebyshev, "Cos2x")
x_value = np.linspace(-np.pi, 3*np.pi, num=500)

error1 = compute_error(x_value, coefs, x_known)
error2 = compute_error(x_value, coefs2, x_chebyshev)
plot_errors(x_value, error1, error2)

