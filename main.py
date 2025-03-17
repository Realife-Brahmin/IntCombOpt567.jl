#%%
# Ensure amplpy is installed (pip install amplpy)
from amplpy import AMPL

# Initialize AMPL environment
ampl = AMPL()

# Define the AMPL model
ampl.eval(
    """
var x {1..4} binary;
# var x {1..4} >= 0;
# var x {1..4} integer >= 0;
maximize z: 17*x[1] + 10*x[2] + 25*x[3] + 17*x[4];

subject to S0:
    5*x[1] + 3*x[2] + 8*x[3] + 7*x[4] <= 12;
"""
)

# Solve the optimization problem
ampl.solve(solver="gurobi")

#%%
# Get results
x_vals_list = ampl.get_variable("x").get_values().to_list()
x_vals = {str(i): val for i, val in x_vals_list}
z_val = ampl.get_objective("z").value()

# Display results
print("Optimal solution:")
for i in range(1, 5):
    print(f"x[{i}] = {x_vals[str(i)]}")

print(f"Optimal objective (z): {z_val}")

# %%
