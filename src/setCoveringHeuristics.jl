module setCoveringHeuristics

export txt2mats

using SparseArrays

function txt2mats(filepath::String)
    # Initialize arrays to store row and column indices for A and A_T
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Read the file and parse the rows
    for line in readlines(filepath)
        i, j = parse.(Int, split(line))
        push!(rows_A, i)
        push!(cols_A, j)
        push!(rows_AT, j)  # Reverse for A_T
        push!(cols_AT, i)
    end

    # Determine the dimensions dynamically
    max_i = maximum(rows_A)
    max_j = maximum(cols_A)

    # Create sparse matrices A and A_T
    A = sparse(rows_A, cols_A, ones(Int, length(rows_A)), max_i, max_j)
    A_T = sparse(rows_AT, cols_AT, ones(Int, length(rows_AT)), max_j, max_i)

    matrixDict = Dict(
        :A => A,
        :A_T => A_T
    )

    return matrixDict
end

end # module setCoveringHeuristics