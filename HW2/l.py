import numpy as np
import matplotlib.pyplot as plt

# Hypothetical error arrays for each root
error_arrays = {
    'root1': np.array([1.16676417e-01, 1.16676417e-01,4.08664478e-02,1.51250129e-02,
 5.71517626e-03, 2.17664503e-03, 8.31481989e-04 ,3.17993267e-04,
 1.21667401e-04 ,4.65590078e-05 ,1.78180923e-05 ,6.81913690e-06,
 2.60976689e-06, 9.98793237e-07]),
    'root2': np.array([6.95139102e-03, 6.95139102e-03, 7.14073829e-04 ,7.87391789e-05,
 8.62146349e-06 ,9.44733640e-07])
}

# Process each errors array
for root, errors in error_arrays.items():
    x = np.log(errors[:-1])
    y = np.log(errors[1:])
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    slopes = dy / dx

    print(f"Slopes for {root} = ", slopes)

    # Plotting
    plt.figure()
    plt.plot(x, y, "bo-", label=f'Log-Log Plot for {root}')
    plt.xlabel("log(e_i)")
    plt.ylabel("log(e_{i+1})")
    plt.title(f"Log-Log Plot of Convergence Errors for {root}")
    plt.grid(True)
    plt.legend()
    # Save to a file
    plt.savefig(f"./LogErrorsPlot_{root}.png")

# Uncomment to display plots (if running interactively)
plt.show()