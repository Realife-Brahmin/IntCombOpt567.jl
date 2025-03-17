# %%
import os
from amplpy import AMPL
ampl = AMPL()

problem = "original"
problem = "S0"
problem = "S1"
problem = "S2"
problem = "S3"
problem = "S4"
problem = "S5"
problem = "S6"
problem = "S7"
# problem = "S8"

# Load AMPL model from external file
model_file_path = os.path.join("rawData", "branchAndBoundExample", "model-"+problem+".mod.txt")
ampl.read(model_file_path)

# Solve the optimization problem
ampl.solve(solver="gurobi")

# Check if the problem was solved
assert ampl.solve_result == "solved"

# Get results
x_vals_list = ampl.get_variable("x").get_values().to_list()
x_vals = {str(i): val for i, val in x_vals_list}
z_val = ampl.get_objective("z").value()

# Display results
print("Optimal solution:")
for i in x_vals.keys():
    print(f"x[{i}] = {x_vals[str(i)]}")

print(f"Optimal objective (z): {z_val}")

# %%
