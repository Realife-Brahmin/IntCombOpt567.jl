module setCoveringHeuristics

export txt2mats

using SparseArrays

function txt2mats(filepath::String)
    # Initialize arrays to store row and column indices for A and A_T
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    degPole = Dict{Int,Int}()  # Degree of each pole (column j)
    degMet = Dict{Int,Int}()   # Degree of each meter (row i)

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
    A_T = sparse(rows_AT, cols_AT, ones(Int, length(rows_AT)), max_j, max_i)

    matrixDict = Dict(
        :A => A,
        :A_T => A_T,
        :degPole => degPole,
        :degMet => degMet
    )

    return matrixDict
end

end # module setCoveringHeuristics