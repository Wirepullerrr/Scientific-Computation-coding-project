import numpy as np
import pandas as pd

# Adjust display settings
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.precision', 15)

# Values of x
x_values = [10**-i for i in range(1, 15)]

# Calculations for Expression 1: Original and Corrected Alternative
original_expr1 = [(1 - 1/np.cos(x)) / np.tan(x)**2 for x in x_values]
corrected_alternative_expr1 = [-1 / (1 + 1/np.cos(x)) for x in x_values]

# Table for Expression 1
results_expr1 = pd.DataFrame({
    'x': x_values,
    'Original Expr1': original_expr1,
    'Alternative Expr1': corrected_alternative_expr1,
})
print("Results for Expression 1:")
print(results_expr1)

