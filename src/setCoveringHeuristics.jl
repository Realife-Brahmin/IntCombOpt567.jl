module setCoveringHeuristics

export txt2graph

using SparseArrays

function txt2graph(filepath::String)
    # Initialize arrays to store row and column indices for A and A_T
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    degPole = Dict{Int, Int}()  # Degree of each pole (column j)
    degMet = Dict{Int, Int}()   # Degree of each meter (row i)

    # Read the file and parse the rows
    for line in readlines(filepath)
        i, j = parse.(Int, split(line))
        push!(rows_A, i)
        push!(cols_A, j)
        push!(rows_AT, j)  # Reverse for A_T
        push!(cols_AT, i)

        # Update degree counts
        degMet[i] = get(degMet, i, 0) + 1
        degPole[j] = get(degPole, j, 0) + 1
    end

    # Determine the dimensions dynamically
    max_i = maximum(rows_A)
    max_j = maximum(cols_A)

    # Create sparse matrices A and A_T
    A = sparse(rows_A, cols_A, ones(Int, length(rows_A)), max_i, max_j)
    A0 = deepcopy(A)
    A_T = sparse(rows_AT, cols_AT, ones(Int, length(rows_AT)), max_j, max_i)
    A_T0 = deepcopy(A_T)

    # Initialize additional data structures
    P = sort(collect(1:max_j))  # 'Set' of all poles
    M = sort(collect(1:max_i))  # 'Set' of all meters
    Pprime = Set{Int}()  # Set of selected poles (initially empty)
    Mprime = Set{Int}()  # Set of covered meters (initially empty)
    Acov = sparse(Int[], Int[], Int[], max_i, max_j)  # Initially empty sparse matrix with known dimensions

    graphState = Dict(
        :A => A,
        :A0 => A0,
        :A_T => A_T,
        :A_T0 => A_T0,
        :Acov => Acov,
        :degPole => degPole,
        :degMet => degMet,
        :M => M,
        :Mprime => Mprime,
        :P => P,
        :Pprime => Pprime,
        :poles_used => 0,
        :meters_covered => 0,
    )

    return graphState
end

end # module setCoveringHeuristics