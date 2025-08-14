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
function read_htp_costs(path::AbstractString; cutoff::Int=10_000, drop_self_loops::Bool=false)
    rows = [parse.(Int, split(strip(line))) for line in eachline(path) if !isempty(strip(line))]
    @assert !isempty(rows) "Empty file: $path"
    n = length(rows)
    @assert all(length(r) == n for r in rows) "Matrix must be square (got row lengths $(length.(rows)))"

    C = fill(Inf, n, n)              # Float64 with Inf for absent arcs
    A = Tuple{Int,Int}[]
    kept_diag = 0
    for i in 1:n, j in 1:n
        cij = rows[i][j]
        if cij < cutoff
            if drop_self_loops && i == j
                continue
            end
            C[i, j] = float(cij)
            push!(A, (i, j))
            kept_diag += (i == j) ? 1 : 0
        end
    end
    N = collect(1:n)
    println("read_htp_costs: n=$n, |A|=$(length(A)), diagonals kept=$(kept_diag > 0)")
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
function solve_htp(path::AbstractString;
    L_list=[200, 500, 800, 1000],
    root::Int=1, mipgap=1e-4, timelimit=0.0,
    case_name::AbstractString="case")   # <-- added
    N, C, A = read_htp_costs(path)
    println("HTP instance loaded: n = $(length(N)), arcs = $(length(A))")

    # NEW: create processedData/<case_name>/ if needed
    outdir = joinpath("processedData", case_name)   # <-- added
    try
        mkpath(outdir)                                      # <-- added
    catch e
        @warn "Could not create output directory $outdir" exception = (e, catch_backtrace())
    end

    for L in L_list
        println("\n=== Solving for L = $L ===")
        model, x, y, f = build_htp_model(N, C, A; L=L, root=root, mipgap=mipgap, timelimit=timelimit)
        optimize!(model)

        term = termination_status(model)
        primal = primal_status(model)
        has_sol = term in (MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.LOCALLY_SOLVED) && primal != MOI.NO_SOLUTION
        @printf("Termination: %-14s  Primal: %-12s\n", string(term), string(primal))
        if !has_sol
            println("No feasible/available solution for L=$L. (Skipping objective/report.)")
            continue
        end

        obj = objective_value(model)
        used_nodes = [i for i in N if value(y[i]) > 0.5]
        num_used_nodes = length(used_nodes)
        total_cost = sum(C[i, j] * value(x[(i, j)]) for (i, j) in A)
        total_arcs = sum(value(x[(i, j)]) for (i, j) in A)

        # NEW: sanity check that model satisfied sum(c*x) ≥ L
        viol = max(0.0, L - obj)
        if viol > 1e-6
            @warn "Model objective is below L by $(viol). This should not happen; printing details."
        end

        @printf("Obj(cost): %.2f   Selected arcs: %.0f   Total cost: %.2f   Visited nodes: %d %s\n",
            obj, total_arcs, total_cost, num_used_nodes, string(used_nodes))

        # NEW: write tour under processedData/<case_name>/
        fname = joinpath(outdir, "tour_L$(L).txt")          # <-- changed
        try
            dump_solution_tour(model, x, A, C, L; filename=fname)
            println("Wrote tour to $fname")
        catch e
            @warn "Could not reconstruct/write tour for L=$L" exception = (e, catch_backtrace())
        end

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
        @printf("Solve time: %s   B&B nodes: %s\n",
            runtime === missing ? "n/a" : @sprintf("%.2f sec", runtime),
            bbnodes === missing ? "n/a" : string(bbnodes))
    end
end



using Printf

# Build an Eulerian circuit from selected arcs (directed)
function eulerian_circuit_from_selected(arcs::Vector{Tuple{Int,Int}}; start::Union{Nothing,Int}=nothing)
    isempty(arcs) && error("No arcs provided to build a tour.")

    # adjacency list with multiplicity
    adj = Dict{Int,Vector{Int}}()
    for (u, v) in arcs
        push!(get!(adj, u, Int[]), v)
        get!(adj, v, Int[])  # ensure key exists
    end

    # choose a start that has outgoing degree
    s = isnothing(start) ? arcs[1][1] : start
    if !haskey(adj, s) || isempty(adj[s])
        # pick any node with outgoing edge
        s = first(k for (k, vs) in adj if !isempty(vs))
    end

    stack = [s]
    path = Int[]
    edges_used = 0
    total_edges = length(arcs)

    while !isempty(stack)
        v = stack[end]
        if haskey(adj, v) && !isempty(adj[v])
            w = pop!(adj[v])   # consume one outgoing arc v->w
            push!(stack, w)
            edges_used += 1
        else
            push!(path, pop!(stack))
        end
    end

    # Path lists vertices; successive pairs are edges of the circuit
    if length(path) < 2
        error("Failed to construct a circuit (path too short).")
    end
    tour_edges = [(path[k], path[k+1]) for k in 1:length(path)-1]

    # In a proper Eulerian circuit, we consume all selected arcs exactly once
    if length(tour_edges) != total_edges
        error("Could not build a single Eulerian circuit; got $(length(tour_edges)) of $total_edges edges. Graph may be disconnected or constraints inconsistent.")
    end
    return tour_edges, path
end

# Write the tour to a text file with i, j, c_ij, cumulative cost
function write_tour_txt(filename::AbstractString, tour_edges::Vector{Tuple{Int,Int}}, C::AbstractMatrix, L::Int)
    open(filename, "w") do io
        cum = 0.0
        @printf(io, "%-6s %-6s %-8s %-12s\n", "i", "j", "c_ij", "cum_cij")
        for (i, j) in tour_edges
            cij = C[i, j]
            if !isfinite(cij)
                error("Arc ($i,$j) is absent in the cost matrix (c_ij = Inf).")
            end
            cum += cij
            @printf(io, "%-6d %-6d %-8.0f %-12.0f\n", i, j, cij, cum)
        end
        first_i = tour_edges[1][1]
        last_j = tour_edges[end][2]
        @printf(io, "\ncloses: %s (start i=%d, end j=%d)\n", last_j == first_i, first_i, last_j)
        @printf(io, "total_cost=%.0f, L=%d, meets_L=%s\n", cum, L, cum >= L)
    end
    return nothing
end

# One-call wrapper: extract selected arcs from (model,x), build circuit, write file
function dump_solution_tour(model, x, A::Vector{Tuple{Int,Int}}, C::AbstractMatrix, L::Int;
    filename::AbstractString="tour.txt", tol::Float64=0.5, start::Union{Nothing,Int}=nothing)
    chosen = [(i, j) for (i, j) in A if value(x[(i, j)]) > tol]
    isempty(chosen) && error("No arcs selected in solution.")
    tour_edges, path = eulerian_circuit_from_selected(chosen; start=start)
    write_tour_txt(filename, tour_edges, C, L)
    return tour_edges, path
end

# ------------------------------
# Example usage (uncomment):
# case = "htp3";L_list = [3, 6, 7];   # expect feasible, feasible, infeasible
# case = "htp_small";L_list = [10, 12, 15, 16, 17];   # expect feasible, feasible, feasible, feasible, infeasible

case = "htp10";L_list=[200,500,800,1000]
# case = "htp20";L_list=[200,500,800,1000]

caseAddr = joinpath("rawData", "project01", case * ".txt") 
solve_htp(caseAddr; L_list=L_list, root=1, mipgap=1e-4, timelimit=60.0, case_name=case)
# ------------------------------
