import numpy as np

def f(x):
    # Cluster sum for k=2
    # P(v) = x
    # Neighbors contribute 2 * (1 - sqrt(x/(1-x)))
    # Total = x + 2 - 2 * sqrt(x/(1-x))
    
    term = np.sqrt(x/(1-x))
    return x + 2 - 2 * term

print("Checking f(x) for x in (1/3, 0.5]")
x_vals = np.linspace(0.334, 0.5, 1000)
y_vals = f(x_vals)

max_y = np.max(y_vals)
max_x = x_vals[np.argmax(y_vals)]

print(f"Max f(x) = {max_y:.5f} at x = {max_x:.5f}")

if max_y <= 1.0 + 1e-9:
    print("PROOF HOLDS: Cluster Sum <= 1 for k=2.")
else:
    print("PROOF FAILS: Cluster Sum > 1 for k=2.")
