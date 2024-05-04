#!/usr/bin/env python
"""
Solve ax^2 + bx + c = 0 for real or complex roots

return 0 in no error, 1 otherwise
"""

import argparse
from math import sqrt, fabs
from cmath import sqrt as csqrt

parser = argparse.ArgumentParser(description="Solve ax^2 + bx + c =0 for real or complex roots")
parser.add_argument("-d","--debug",action="store_true")
DEBUG = parser.parse_args().debug
r1 = None
r2 = None
############################## FUNCTIONS #############################

def quadraticFormula(a,b,c):
    global r1, r2
    discriminant = b*b - 4*a*c
    eps = 1e-20

    # Define a tolerance for considering the discriminant effectively zero
    tolerance = 1e-14

    if -tolerance < discriminant < tolerance:
        # Treat the discriminant as zero
        discriminant = 0

    if DEBUG:
        print("a = %.20f" % a)
        print("b = %.20f" % b)
        print("c = %.20f" % c)
        print("D = %.20f" % discriminant)
    
    if(discriminant < 0 and fabs(discriminant) < eps):
        print("|abs(D)| =  %.6e; Setting D to 0" % fabs(discriminant))
        discriminant = 0
    
    if fabs(a) < eps:
        return 1
    
    if discriminant > 0:
        r1 = (-b + sqrt(discriminant)) / (2*a)
        r2 = (-b - sqrt(discriminant)) / (2*a)
        print("Real and distinct roots")
    elif discriminant == 0:
        r1 = r2 = -b / (2*a)
        print("Real and equal roots")
    else:
        r1 = (-b + csqrt(discriminant)) / (2*a)
        r2 = (-b - csqrt(discriminant)) / (2*a)
        print("Complex roots")
    return 0

############################## MAIN #############################

print("Solve ax^2 + bx + c = 0 for real or complex roots.")
a = float(input("Enter a: "))
b = float(input("Enter b: "))
c = float(input("Enter c: "))

error = quadraticFormula(a,b,c)


if error:
    print("ERROR:Invalid inputs for a quadratic equation, such as a=0")
else:
    print("Roots are {} and {}".format(r1, r2))