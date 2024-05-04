import numpy as np
import pandas as pd

# Adjust display settings for more columns, rows, and precision
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.precision', 15)  # Display up to 15 decimal places

# Values of x
x_values = [10**-i for i in range(1, 15)]


# Calculations for Expression 2: Original and Alternative
original_expr2 = [(1 - (1 - x)**3) / x for x in x_values]
alternative_expr2 = [3 - 3*x + x**2 for x in x_values]

# Table for Expression 2
results_expr2 = pd.DataFrame({
    'x': x_values,
    'Original Expr2': original_expr2,
    'Alternative Expr2': alternative_expr2,
})

# Print results for Expression 2
print("Results for Expression 2:")
print(results_expr2)
