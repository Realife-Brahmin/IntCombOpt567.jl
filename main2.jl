using JuMP
using Gurobi
using LinearAlgebra
using Printf

"""
    read_htp_costs(path::AbstractString; cutoff=10_000)

Read the HTP cost matrix from `path`. Returns:
- N::Vector{Int}   : node indices 1..n
- C::Matrix{Int}   : costs (cij), with `Inf` where arc is absent
- A::Vector{Tuple{Int,Int}} : list of allowed directed arcs (i,j)
"""
function read_htp_costs(path::AbstractString; cutoff=10_000)
    raw = readlines(path)
    data = [parse.(Int, split(strip(line))) for line in raw if !isempty(strip(line))]
    Cint = reduce(hcat, data)'  # rows = i, cols = j
    n = size(Cint, 1)
    N = collect(1:n)
    C = fill(Inf, n, n)
    A = Tuple{Int,Int}[]
    for i in N, j in N
        cij = Cint[i, j]
        if cij < cutoff
            C[i, j] = cij
            push!(A, (i, j))
        end
    end
    return N, C, A
end

"""
    build_htp_model(N, C, A; L::Int, root::Int=1, mipgap=1e-4, timelimit=0.0)

Build a JuMP model for the Hiker's Tour Problem:

- Decision vars:
    x[i,j] ∈ {0,1}  (use arc)
    y[i]   ∈ {0,1}  (node used)
    f[i,j] ≥ 0      (connectivity flow, capacity-limited by x)
- Constraints:
    1) Min total length ≥ L
    2) In-degree == Out-degree (Eulerian balance) at every node
    3) Node-arc linking: if no incident arcs, y[i]=0; if arcs used, y[i]=1
    4) Single connected component among used nodes (single-commodity flow)
    5) Non-empty: sum(y) ≥ 1
- Objective: minimize total cost ∑ c_ij x_ij subject to (1)–(5)

Returns (model, x, y, f).
"""
function build_htp_model(N, C, A; L::Int, root::Int=1, mipgap=1e-4, timelimit=0.0)
    n = length(N)
    @assert 1 ≤ root ≤ n "root must be a valid node index"

    env = Gurobi.Env()
    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    # Optional controls
    set_optimizer_attribute(model, "MIPGap", mipgap)
    if timelimit > 0
        set_optimizer_attribute(model, "TimeLimit", timelimit)
    end

    # Convenience index sets per node
    out_arcs = Dict(i => [(i,j) for (i2,j) in A if i2 == i] for i in N)
    in_arcs  = Dict(i => [(j,i) for (j,i2) in A if i2 == i] for i in N)

    # Decision variables
    @variable(model, x[A], Bin)                     # select arc
    @variable(model, y[N], Bin)                     # node used
    @variable(model, f[A] >= 0.0)                   # flow for connectivity

    # Objective: minimize total cost while meeting L
    @objective(model, Min, sum(C[i,j] * x[(i,j)] for (i,j) in A))

    # 1) Minimum tour length
    @constraint(model, sum(C[i,j] * x[(i,j)] for (i,j) in A) >= L)

    # 2) Eulerian in/out balance at each node
    @constraint(model, [i in N], sum(x[a] for a in out_arcs[i]) - sum(x[a] for a in in_arcs[i]) == 0)

    # 3) Node-arc linking (if node unused, no incident arcs are used; if any arc used, node is used)
    # Use big-M equal to node's max possible degree
    for i in N
        M_i_out = length(out_arcs[i])
        M_i_in  = length(in_arcs[i])
        @constraint(model, sum(x[a] for a in out_arcs[i]) <= M_i_out * y[i])
        @constraint(model, sum(x[a] for a in in_arcs[i])  <= M_i_in  * y[i])
        # If either indeg or outdeg is positive, y[i] must be 1
        @constraint(model, y[i] <= (sum(x[a] for a in out_arcs[i]) + sum(x[a] for a in in_arcs[i])))
    end

    # 4) Single connected component among used nodes via single-commodity flow
    #    Root supplies ∑_{i≠root} y[i] units. Each used non-root node consumes 1 unit.
    #    Capacities: f[i,j] ≤ U * x[i,j], with U = n-1 (enough to reach all nodes).
    U = n - 1
    @constraint(model, [a in A], f[a] <= U * x[a])

    @constraint(model, [i in N; i != root],
        sum(f[(j, i)] for (j, i) in in_arcs[i]) -
        sum(f[(i, j)] for (i, j) in out_arcs[i]) == y[i])


    @constraint(model,
        sum(f[(root,j)] for (root,j) in out_arcs[root]) - sum(f[(i,root)] for (i,root) in in_arcs[root]) ==
        sum(y[i] for i in N if i != root))

    # 5) Non-empty tour (can also ensure at least one arc chosen)
    @constraint(model, sum(y) >= 1)

    return model, x, y, f
end

"""
    solve_htp(path::AbstractString; L_list=[200,500,800,1000], root=1, mipgap=1e-4, timelimit=0.0)

Reads the instance, solves for each L, and prints a small report.
"""
function solve_htp(path::AbstractString; L_list=[200, 500, 800, 1000], root=1, mipgap=1e-4, timelimit=0.0)
    N, C, A = read_htp_costs(path)
    println("HTP instance loaded: n = $(length(N)), arcs = $(length(A))")

    for L in L_list
        println("\n=== Solving for L = $L ===")
        model, x, y, f = build_htp_model(N, C, A; L=L, root=root, mipgap=mipgap, timelimit=timelimit)
        optimize!(model)

        term = termination_status(model)
        primal = primal_status(model)
        obj   = objective_value(model)
        used_nodes = [i for i in N if value(y[i]) > 0.5]
        num_used_nodes = length(used_nodes)
        total_cost = sum(C[i,j] * value(x[(i,j)]) for (i,j) in A)
        total_arcs = sum(value(x[(i,j)]) for (i,j) in A)

        # Gurobi stats
        bbnodes = try
            MOI.get(model, MOI.NodeCount())
        catch
            missing
        end
        runtime = try
            MOI.get(model, MOI.SolveTime())
        catch
            missing
        end

        @printf("Status: %-12s  Primal: %-10s  Obj(cost): %.0f\n", string(term), string(primal), obj)
        @printf("Visited nodes: %d  (indices: %s)\n", num_used_nodes, string(used_nodes))
        @printf("Selected arcs: %.0f  Total cost (from x): %.0f  (constraint L=%d)\n", total_arcs, total_cost, L)
        @printf("Solve time: %s sec   B&B nodes: %s\n",
                isnothing(runtime) || runtime === missing ? "n/a" : @sprintf("%.2f", runtime),
                bbnodes === missing ? "n/a" : string(bbnodes))
    end
end

# ------------------------------
# Example usage (uncomment):
case = "htp_small"
case = "htp10"
case = "htp20"
solve_htp("rawData/project01/"*case*".txt"; L_list=[200,500,800,1000], root=1, mipgap=1e-4, timelimit=60.0)
# ------------------------------
