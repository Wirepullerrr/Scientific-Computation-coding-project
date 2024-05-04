# Newton Divided Difference Interpolation

import numpy as np
import matplotlib.pyplot as pyp
import argparse
# this just makes it so that if you do "python3 newtonDD.py --test", it will run test environment
# where numbers are input through the command line. You can modify how this is carried out by
# changing the testNewton() function
'''parser = argparse.ArgumentParser(description="Provide Newton's Divided Difference Interpolation")
parser.add_argument("-t","--test",action="store_true")
TEST = parser.parse_args().test'''
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

# x,y are 1d- arrays
def newtonDD(x,y):
    n = len(x)
    if (len(y) != n): exit(1)
    coefs = newtonDDsetup(x, y)
    printPolynomial(coefs, x)
    x_values = np.linspace(np.min(x), np.max(x), 500)
    y_values = [newtonEval(x, y, coefs) for x in x_values]

    estimate_2010 = newtonEval(2010, coefs, x)
    print()
    print("Estimate 2010 value = " + str(estimate_2010))
    '''if TEST: print("x =",x)
    if TEST: print("y =",y)
    
    if TEST: 
        print("The coefs are: ", end=" ")
        for i in range(n):
            print(" %g"%(coefs[i]),end=" ")
        print()'''

    m = 10*n
    minx = min(x)
    maxx = max(x)
    t = np.arange(minx,maxx+1,(maxx-minx)/m)
    val = newtonEval(t,coefs,x)

    # plot
    pyp.plot(t,val)
    for i in range(n):
        pyp.plot(x[i],y[i],'k*')
    pyp.xlabel("x")
    pyp.ylabel("y")
    pyp.xlim(1994,2010)
    pyp.ylim(-400, 100)
    pyp.title("Newton DD interpolation")
    pyp.legend(["Newton DD interpolant","Data points"],loc="best")
    pyp.show()



# Estimate for 2010 using the interpolating polynomial


def printPolynomial(coefs, x):
    n = len(coefs)

    print("Interpolating Polynomial:")
    print("P(x) =")

    for i in range(n - 1, -1, -1):  # Start from the last coefficient
        term = f"{coefs[i]:+.6f}"  # Include sign in format
        for j in range(i):
            term += f" * (x - {x[j]:.0f})"
        if i > 0:
            print(f"{term}")
        else:
            print(term)

'''def printPolynomial(coefs, x):
    n = len(coefs)
    terms = []

    # Construct the polynomial as a string
    for i in range(n):
        term = f"{coefs[i]:.6f}"
        for j in range(i):
            term += f"*(x - {x[j]:.3f})"
        terms.append(term)

    # Combine terms into a polynomial string
    polynomial = " + ".join(terms)

    print("Interpolating Polynomial:")
    print(f"P(x) = {polynomial}")'''


'''def testNewton():
    n = int(input("Enter n: "))

    # Get data from command line
    ans = input("Enter data by points? [y/n] ")
    if (ans[0] == 'y' or ans[0] == 'Y'):
        data = input(f"Enter {n} data points (x1 y1 x2 y2 ... xn yn): ").split(" ")
        # other ways to do this; just for testing without entering a file.
        # pull all the x-values and convert to floats: even indices 0, 2, 4, ...
        x = list(map(float, np.array(data)[np.arange(0,len(data),2)]))
        # pull all the y-values and convert to floats: odd indices 1, 3, 5, ...
        y = list(map(float, np.array(data)[np.arange(1,len(data),2)]))
    else:
        x_data = input(f"Enter {n} distinct x values (x1 x2 ... xn): ").split(" ")
        x = list(map(float,x_data))
        y_data = input(f"Enter {n} distinct y values (y1 y2 ... yn): ").split(" ")
        y = list(map(float,y_data))
    
    newtonDD(x,y)'''
    




############################## MAIN #############################
years = np.array([1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003])
production = np.array([67.052, 68.008, 69.803, 72.024, 73.400, 72.063, 74.669, 74.487, 74.065, 76.777])
newtonDD(years, production)



