from sympy import symbols, sqrt, diff, solve, simplify

x = symbols('x', real=True)
# Constrain x to (1/3, 1/2]
# f(x) = x + 2 - 2 * sqrt(x / (1-x))
f = x + 2 - 2 * sqrt(x / (1 - x))

print("Function f(x):", f)

# Derivative
df = diff(f, x)
print("Derivative f'(x):", simplify(df))

# Evaluate at boundaries
val_1_3 = f.subs(x, 1/3).evalf()
val_1_2 = f.subs(x, 1/2).evalf()

print(f"f(1/3) = {val_1_3}")
print(f"f(1/2) = {val_1_2}")

# Find critical points where derivative is 0
crit = solve(df, x)
print("Critical points:", crit)
