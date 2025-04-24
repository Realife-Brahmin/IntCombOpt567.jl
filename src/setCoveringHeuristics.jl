module setCoveringHeuristics

export
    chooseNextPole,
    selectPole,
    solveSetCoveringProblem,
    txt2graph

using Parameters
using SparseArrays

include("./helperFunctions.jl")
import .helperFunctions as HF

function txt2graph(filepath::String)
    # Initialize arrays to store row and column indices for A and A_T
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    degPoleUnused = Dict{Int, Int}()  # Degree of each pole (column j)
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
        degPoleUnused[j] = get(degPoleUnused, j, 0) + 1
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
    p = length(P)
    M = sort(collect(1:max_i))  # 'Set' of all meters
    m = length(M)
    Pprime = Set{Int}()  # Set of selected poles (initially empty)
    Mprime = Set{Int}()  # Set of covered meters (initially empty)
    Acov = sparse(Int[], Int[], Int[], max_i, max_j)  # Initially empty sparse matrix with known dimensions

    graphState = Dict(
        :A => A,
        :A0 => A0,
        :A_T => A_T,
        :A_T0 => A_T0,
        :Acov => Acov,
        :cleanupDoneLastIter => false,
        :cleanupUsefulLastIter => false,
        :degPoleUnused => degPoleUnused,
        :degMet => degMet,
        :M => M,
        :m => m,
        :Mprime => Mprime,
        :P => P,
        :p => p,
        :Pprime => Pprime,
        :poles_used => 0,
        :meters_covered => 0,
    )

    return graphState
end

function chooseNextPole(graphState)
    @unpack degPoleUnused = graphState
    if isempty(degPoleUnused)
        error("No unused poles available.")
    end
    j_candidate = HF.argmax_smallestkey(degPoleUnused)[1]  # Pole with the maximum degree
    return j_candidate
end

function selectPole(graphState, j;
    verbose = false)
    @unpack A, A_T, degPoleUnused, degMet, Mprime, Pprime, Acov = graphState

    # Update the set of covered meters (Mprime) and selected poles (Pprime)

    HF.myprintln(verbose, "Pole $j selected")
    meters_covered_by_j = findall(A[:, j] .== 1)  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")
    Mprime = union(Mprime, meters_covered_by_j)  # Add these meters to Mprime
    Pprime = union(Pprime, j)  # Add pole j to Pprime

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j
        degMet[i] -= 1
    end

    # Remove pole j from degPoleUnused
    delete!(degPoleUnused, j)

    # Update the sparse matrix Acov to reflect the meters covered by pole j
    for i in meters_covered_by_j
        Acov[i, j] = 1  # Set Acov[i, j] = 1
    end

    # Remove pole j from A (set A[i, j] = 0 for all i)
    for i in meters_covered_by_j
        A[i, j] = 0
    end

    poles_used = length(Pprime)  # Update the number of poles used
    # Update the graph state
    @pack! graphState = Acov, Mprime, Pprime, poles_used, A, degPoleUnused, degMet
    return graphState
end

function solveSetCoveringProblem(graphState;
    verbose::Bool = false)
    @unpack m, Mprime = graphState # Mprime initially is empty
    
    k = 0
    shouldStop = false
    while !shouldStop # While there are still uncovered meters
        k += 1
        HF.myprintln(verbose, "Iteration $(k): Currently covered meters: $(Mprime)")
        j = chooseNextPole(graphState)  # Choose the next pole
        graphState = selectPole(graphState, j)  # Select the pole and update the graph state
        @unpack Mprime = graphState

        shouldStop = checkForStoppingCriteria(graphState)
    end

    return graphState
end

function checkForStoppingCriteria(graphState;
    verbose::Bool = false)
    @unpack m, Mprime = graphState
    shouldStop = false
    allMetersCovered = false  
    if  length(Mprime) == m  # All meters are covered
        HF.myprintln(true, "All meters are covered.")
        allMetersCovered = true
    end

    if allMetersCovered
        HF.myprintln(true, "Stopping criteria met: All meters are covered.")
        shouldStop = true
    end

    return shouldStop

end

end # module setCoveringHeuristics