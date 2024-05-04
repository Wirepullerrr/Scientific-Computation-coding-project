import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def f(x):
    return (1 + x) ** (-1)

def f_prime(x):
    return -(1 + x) ** (-2)


x = 1
exact_derivative = f_prime(x)

h = np.logspace(-1, -12, num=12, base=10)

# Calculate derivatives using the three-point centered-difference formula
approx_derivatives = (f(x + h) - f(x - h)) / (2 * h)

# Calculate errors
errors = np.abs(approx_derivatives - exact_derivative)

data = {'h': h, 'Approximate f\'(x)': approx_derivatives, 'Error': errors}
data_frame = pd.DataFrame(data)
print(data_frame)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.loglog(h, errors, marker='o')
plt.xlabel('h')
plt.ylabel('Error')
plt.title('Log-Log Plot of Error vs. h for the Three-Point Centered-Difference Formula')
plt.grid(True, which="both", ls="--")
plt.show()


