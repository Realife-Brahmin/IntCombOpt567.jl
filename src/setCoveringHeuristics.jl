module setCoveringHeuristics

export
    addPole!,
    checkForRedundantPole,
    checkForStoppingCriteria,
    cleanupGraph!,
    chooseNextPole,
    initializeGraph,
    solveSetCoveringProblem!,
    removePole!

using Parameters
using SparseArrays

include("./helperFunctions.jl")
import .helperFunctions as HF

function initializeGraph(filepath::String;
    maxiter::Int = 100000)
    # Initialize arrays to store row and column indices for A and A_T
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    degPoleUnused = Dict{Int, Int}()  # Degree of each pole (column j)
    degMetUnusedPoles = Dict{Int, Int}()   # Degree of each meter (row i) wrt poles NOT in Pprime
    degMetUsedPoles = Dict{Int, Int}()  # Degree of each meter (row i) wrt poles part of Pprime

    # Read the file and parse the rows
    for line in readlines(filepath)
        i, j = parse.(Int, split(line))
        push!(rows_A, i)
        push!(cols_A, j)
        push!(rows_AT, j)  # Reverse for A_T
        push!(cols_AT, i)

        # Update degree counts
        degMetUnusedPoles[i] = get(degMetUnusedPoles, i, 0) + 1
        degMetUsedPoles[i] = 0
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
        :cleanupRepeats => 10,
        :cleanupUsefulLastIter => false,
        :degPoleUnused => degPoleUnused,
        :degMetUsedPoles => degMetUsedPoles,
        :degMetUnusedPoles => degMetUnusedPoles,
        :k => 0,
        :M => M,
        :m => m,
        :maxiter => maxiter,
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

function addPole!(graphState, j;
    verbose = false)
    @unpack A, A_T, A_T0, degPoleUnused, degMetUsedPoles, degMetUnusedPoles, Mprime, Pprime, Acov = graphState

    HF.myprintln(verbose, "Pole $j to be added")
    meters_covered_by_j = findall(A_T0[j, :] .== 1)  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")
    Mprime = union(Mprime, meters_covered_by_j)  # Add these meters to Mprime
    Pprime = union(Pprime, j)  # Add pole j to Pprime

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j
        degMetUnusedPoles[i] -= 1
        degMetUsedPoles[i] += 1
    end

    # Remove pole j from degPoleUnused
    delete!(degPoleUnused, j)

    # Update the sparse matrix Acov to reflect the meters covered by pole j
    for i in meters_covered_by_j
        Acov[i, j] = 1 
    end

    # Remove pole j from A (set A[i, j] = 0 for all i)
    for i in meters_covered_by_j
        A[i, j] = 0
    end
    dropzeros!(A)

    poles_used = length(Pprime)
    meters_covered = length(Mprime) 
    # Update the graph state
    @pack! graphState = Acov, Mprime, Pprime, poles_used, A, degPoleUnused, degMetUsedPoles, degMetUnusedPoles, meters_covered
    return graphState
end

function removePole!(graphState, j;
    verbose::Bool = false)
    @unpack A, A_T, A_T0, degPoleUnused, degMetUsedPoles, degMetUnusedPoles, Mprime, Pprime, Acov = graphState

    if j ∉ Pprime
        error("Attempting to remove a pole that is not in P′.")
        return
    end

    Pprime = setdiff(Pprime, j)  # Remove pole j from Pprime

    HF.myprintln(verbose, "Pole $j to be removed")
    meters_covered_by_j = findall(A_T0[j, :] .== 1)  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j
        degMetUnusedPoles[i] += 1  # Update degrees for meters covered by pole j
        degMetUsedPoles[i] -= 1
        if degMetUsedPoles[i] == 0
            Mprime = setdiff(Mprime, i)  # Remove meter i from Mprime if it is no longer covered by any pole
        end
    end

    # Update the sparse matrix Acov to reflect the meters no longer covered by pole j
    for i in meters_covered_by_j
        Acov[i, j] = 0 
    end
    dropzeros!(Acov)
    # Add back pole j to A (set A[i, j] = 1 for all i)
    for i in meters_covered_by_j
        A[i, j] = 1
    end

    poles_used = length(Pprime)
    meters_covered = length(Mprime)
    # Update the graph state
    @pack! graphState = Acov, Mprime, Pprime, poles_used, A, degPoleUnused, degMetUsedPoles, degMetUnusedPoles, meters_covered
    return graphState

end

function solveSetCoveringProblem!(graphState;
    verbose::Bool = false)
    @unpack cleanupRepeats, m = graphState
    shouldStop = false
    while !shouldStop # While there are still uncovered meters
        @unpack k, Mprime = graphState # Starts at 0
        k += 1  # Increment the iteration count
        HF.myprintln(verbose, "Iteration $(k): Currently covered meters: $(Mprime)")
        j = chooseNextPole(graphState)  # Choose the next pole
        addPole!(graphState, j)  # Select the pole and update the graph state

        @unpack meters_covered = graphState
        @pack! graphState = k # k-th iteration completed, so saving it

        cleanupAttempt = false
        if meters_covered == m
            HF.myprintln(verbose, "Iteration $(k): Attempting cleanup since all meters are covered")
            cleanupAttempt = true
        elseif k % cleanupRepeats == 0
            HF.myprintln(verbose, "Iteration $(k): Attempting periodic cleanup procedure")
            cleanupAttempt = true
        end
            
        if cleanupAttempt
            cleanupGraph!(graphState; verbose=verbose)
        end

        shouldStop = checkForStoppingCriteria(graphState)
    end

    return graphState
end

function checkForStoppingCriteria(graphState;
    verbose::Bool = false)
    @unpack m, Mprime, k, maxiter = graphState

    if k >= maxiter
        HF.myprintln(true, "Maximum iterations reached!")
        return true
    end

    if length(Mprime) == m
        HF.myprintln(true, "Stopping criterion met: All meters are covered")
        return true
    end

    return false
end

function cleanupGraph!(graphState;
    verbose::Bool = false)

    @unpack poles_used, Pprime = graphState
    cleanupUsefulLastIter = false

    PprimeSorted = sort(collect(Pprime))  # Sort the poles in Pprime

    for j ∈ PprimeSorted
        redundant = checkForRedundantPole(graphState, j; verbose=verbose)
        if redundant
            HF.myprintln(verbose, "Pole $j is redundant and will be removed")
            removePole!(graphState, j; verbose=verbose)
            cleanupUsefulLastIter = true
        end
    end

    @pack! graphState = cleanupUsefulLastIter # Update the graph state with the cleanup result
    return graphState
end

function checkForRedundantPole(graphState, j;
    verbose::Bool = false)

    HF.myprintln(verbose, "Checking if pole $j is redundant")
    @unpack degMetUsedPoles, A_T0 = graphState
    meters_covered_by_j = findall(A_T0[j, :] .== 1)  # Find all meters covered by pole j
    
    for i in meters_covered_by_j
        if degMetUsedPoles[i] == 1  # If meter i is only covered by pole j
            HF.myprintln(verbose, "Pole $j is NOT redundant as it covers meter $i exclusively")
            return false
        end
    end

    return true
end

end # module setCoveringHeuristics