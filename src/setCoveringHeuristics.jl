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
    maxiter::Int = 100000,
    cleanupRepeats::Int = 10,
    scoring_function::String = "greedy")
    # Initialize arrays to store row and column indices for A_m2p_remaining and A_p2m_uncovered
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    degMetUncovered = Dict{Int, Int}()  # Degree of each meter not currently in Mprime
    degPolesRemaining = Dict{Int, Int}()  # Degree of each pole (column j)
    degMetUnusedPoles = Dict{Int, Int}()   # Degree of each meter (row i) wrt poles NOT in Pprime
    degMetUsedPoles = Dict{Int, Int}()  # Degree of each meter (row i) wrt poles part of Pprime

    # Read the file and parse the rows
    for line in readlines(filepath)
        i, j = parse.(Int, split(line))
        push!(rows_A, i)
        push!(cols_A, j)
        push!(rows_AT, j)  # Reverse for A_p2m_uncovered
        push!(cols_AT, i)

        # Update degree counts
        degMetUnusedPoles[i] = get(degMetUnusedPoles, i, 0) + 1
        degMetUncovered[i] = degMetUnusedPoles[i]
        degMetUsedPoles[i] = 0
        degPolesRemaining[j] = get(degPolesRemaining, j, 0) + 1
    end

    degPoleUncoveredMeters = deepcopy
    # Determine the dimensions dynamically
    max_i = maximum(rows_A)
    max_j = maximum(cols_A)

    # Create sparse matrices A_m2p_remaining and A_p2m_uncovered
    A_m2p_remaining = sparse(rows_A, cols_A, ones(Int, length(rows_A)), max_i, max_j)
    A_m2p_reference = deepcopy(A_m2p_remaining)
    A_p2m_uncovered = sparse(rows_AT, cols_AT, ones(Int, length(rows_AT)), max_j, max_i)
    A_p2m_ref = deepcopy(A_p2m_uncovered)

    Aadj_m2p_ref = Dict{Int,Vector{Int}}()
    Aadj_m2p = Dict{Int,Vector{Int}}()
    Aadj_p2m_ref = Dict{Int,Vector{Int}}()
    Aadj_p2m = Dict{Int,Vector{Int}}()

    # Initialize adjacency lists
    for i in 1:max_i
        Aadj_m2p[i] = Vector{Int}()
        Aadj_m2p_ref[i] = Vector{Int}()
    end

    for j in 1:max_j
        Aadj_p2m_ref[j] = Vector{Int}()
        Aadj_p2m[j] = Vector{Int}()
    end

    # Populate adjacency lists
    for (i, j) in zip(rows_A, cols_A)
        push!(Aadj_m2p_ref[i], j)  # Add pole j to the list of poles for meter i
        push!(Aadj_p2m_ref[j], i)  # Add meter i to the list of meters for pole j
    end

    # Initialize additional data structures
    P = sort(collect(1:max_j))  # 'Set' of all poles
    Premaining = deepcopy(P)
    p = length(P)
    M = sort(collect(1:max_i))  # 'Set' of all meters
    m = length(M)
    Pprime = Set{Int}()  # Set of selected poles (initially empty)
    Mprime = Set{Int}()  # Set of covered meters (initially empty)
    A_m2p = sparse(Int[], Int[], Int[], max_i, max_j)  # Initially empty sparse matrix with known dimensions
    A_p2m = sparse(Int[], Int[], Int[], max_j, max_i)  # Initially empty sparse matrix with known dimensions

    graphState = Dict(
        :A_m2p => A_m2p,
        :Aadj_m2p => Aadj_m2p,
        :Aadj_m2p_ref => Aadj_m2p_ref,
        :A_m2p_remaining => A_m2p_remaining,

        :A_p2m => A_p2m,
        :Aadj_p2m => Aadj_p2m,
        :Aadj_p2m_ref => Aadj_p2m_ref,
        :A_p2m_uncovered => A_p2m_uncovered,

        :cleanupDoneLastIter => false,
        :cleanupRepeats => cleanupRepeats,
        :cleanupUsefulLastIter => false,

        :degPolesRemaining => degPolesRemaining,
        :degMetUncovered => degMetUncovered,
        :degMetUsedPoles => degMetUsedPoles,
        :degMetUnusedPoles => degMetUnusedPoles,
        
        :k => 0,
        :M => M,
        :m => m,
        :maxiter => maxiter,
        :Mprime => Mprime,
        :P => P,
        :Premaining => Premaining,
        :p => p,
        :Pprime => Pprime,
        :poles_used => 0,
        :meters_covered => 0,
        :scoring_function => scoring_function,
    )

    return graphState
end

function chooseNextPole(graphState; verbose::Bool = false)

    @unpack degPolesRemaining, scoring_function = graphState
    if isempty(degPolesRemaining)
        error("No unused poles available.")
    end

    if scoring_function == "score1" || scoring_function == "score2"
        @unpack degMetUncovered, Aadj_m2p_ref = graphState
        k = parse(Int, last(scoring_function))  # Extract the number from "score1" or "score2"
        scoreDict = compute_score(degMetUncovered, degPolesRemaining, Aadj_m2p_ref, k=k, verbose=verbose)
    elseif scoring_function == "greedy"
        scoreDict = degPolesRemaining 
    else
        @error("Invalid scoring function: $scoring_function")
    end

    j_candidate = HF.argmax_smallestkey(scoreDict)[1]  # Pole with the maximum degree
    # j_candidate = argmax(scoreDict)[1]  # Pole with the maximum degree

    return j_candidate
end

function addPole!(graphState, j;
    verbose = false)
    @unpack A_m2p_remaining, Aadj_m2p_ref,  Aadj_p2m_ref, degPolesRemaining, degMetUncovered, degMetUsedPoles, degMetUnusedPoles, Mprime, Pprime, Premaining, A_m2p = graphState

    HF.myprintln(verbose, "Pole $j to be added")
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")

    Pprime = union(Pprime, j)  # Add pole j to Pprime
    Premaining = setdiff(Premaining, j)  # Selected pole no longer on market
    delete!(degPolesRemaining, j) # Selected pole no longer on market

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j
        degMetUnusedPoles[i] -= 1 # one less remaining pole which can cover meter i
        degMetUsedPoles[i] += 1 # another pole which now covers meter i
        if !(i ∈ Mprime) # i.e. a previously uncovered meter is being covered by pole j
            delete!(degMetUncovered, i)  # meter i now no longer uncovered

            # All poles covering meter i still available in the market (excluding newly added pole j, already removed from Premaining)
            other_remaining_poles_also_covering_i = intersect(Aadj_m2p_ref[i], Premaining)

            for j_other in other_remaining_poles_also_covering_i
                degPolesRemaining[j_other] -= 1  # Now that a meter is covered without pole j_other's help, its degree is deducted by 1
            end
        end
    end

    Mprime = union(Mprime, meters_covered_by_j)  # Add these meters to Mprime
    # Modifying Mprime only now as we need to check if the meters are already in Mprime


    # Update the sparse matrix A_m2p to reflect the meters covered by pole j
    @unpack A_m2p = graphState
    for i in meters_covered_by_j
        A_m2p[i, j] = 1 
    end
    @pack! graphState = A_m2p

    # Remove pole j from A_m2p_remaining (set A_m2p_remaining[i, j] = 0 for all i)
    @unpack A_m2p_remaining = graphState
    for i in meters_covered_by_j
        A_m2p_remaining[i, j] = 0
    end
    @pack! graphState = A_m2p_remaining

    poles_used = length(Pprime)
    meters_covered = length(Mprime) 
    # Update the graph state
    @pack! graphState = A_m2p, Mprime, Pprime, poles_used, Premaining, degPolesRemaining, degMetUncovered, degMetUsedPoles, degMetUnusedPoles, meters_covered
    return graphState
end

function compute_score(degMetUncovered, degPolesRemaining, Aadj_m2p_ref; verbose::Bool = false, k=1)

    htc_threshold = minimum(values(degMetUncovered))  # Hard-to-cover threshold (minimum degree of uncovered meters)
    HF.myprintln(verbose, "Hard-to-cover threshold: $htc_threshold")
    scoreDict = Dict{Int, Int}()
    
    if k == 1
        for i ∈ keys(degMetUncovered)
            if degMetUncovered[i] == htc_threshold
                # HF.myprintln(true, "Meter $i is a hard to cover meter")
                for j ∈ Aadj_m2p_ref[i]
                    # HF.myprintln(true, "Meter $i is covered by pole $j")
                    scoreDict[j] = degPolesRemaining[j]
                end
            end
        end
    elseif k == 2
        for i ∈ keys(degMetUncovered)
            if degMetUncovered[i] == htc_threshold
                # HF.myprintln(true, "Meter $i is a hard to cover meter")
                for j ∈ Aadj_m2p_ref[i]
                    # HF.myprintln(true, "Meter $i is covered by pole $j")
                    scoreDict[j] = get(scoreDict, j, 0) + degPolesRemaining[j]
                end
            end
        end
    else
        @error("Invalid value for k: $k")
    end

    return scoreDict
end

function removePole!(graphState, j;
    verbose::Bool = false)
    @unpack Aadj_m2p_ref, Aadj_p2m_ref, Premaining, degPolesRemaining, degMetUncovered, degMetUsedPoles, degMetUnusedPoles, Mprime, Pprime  = graphState

    if j ∉ Pprime
        error("Attempting to remove a pole that is not in P′.")
        return
    end

    Pprime = setdiff(Pprime, j)  # Remove pole j from Pprime
    # But j will not be placed back in Premaining (by design)

    HF.myprintln(verbose, "Pole $j to be removed")
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")
    
    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j
        # Aadj_m2p[i] = setdiff(Aadj_m2p[i], j)  # Remove pole j from the adjacency list of meter i
        degMetUnusedPoles[i] += 1  # Update degrees for meters covered by pole j
        degMetUsedPoles[i] -= 1
        if degMetUsedPoles[i] == 0
            @error("Removal of pole $j leaves meter $i uncovered! Why was this selected?")
        end
    end

    # Update the sparse matrix A_m2p to reflect the meters no longer covered by pole j
    @unpack A_m2p = graphState
    for i in meters_covered_by_j
        A_m2p[i, j] = 0 
    end
    @pack! graphState = A_m2p

    
    # Add back pole j to A_m2p_remaining (set A_m2p_remaining[i, j] = 1 for all i)
    @unpack A_m2p_remaining = graphState    
    for i in meters_covered_by_j
        A_m2p_remaining[i, j] = 1
    end
    @pack! graphState = A_m2p_remaining

    poles_used = length(Pprime)
    meters_covered = length(Mprime)
    # Update the graph state
    @pack! graphState = Mprime, Pprime, poles_used, degPolesRemaining, degMetUncovered, degMetUsedPoles, degMetUnusedPoles, meters_covered
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
        HF.myprintln(verbose, "Stopping criterion met: All meters are covered")
        return true
    end

    @unpack poles_used, p = graphState
    if length(poles_used) == p
        HF.myprintln(true, "Stopping criterion met: All poles are used")
        @warn("All poles are used, but not all meters are covered!")
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
    @unpack degMetUsedPoles, Aadj_p2m_ref = graphState
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    
    for i in meters_covered_by_j
        if degMetUsedPoles[i] == 1  # If meter i is only covered by pole j
            HF.myprintln(verbose, "Pole $j is NOT redundant as it covers meter $i exclusively")
            return false
        end
    end

    return true
end

end # module setCoveringHeuristics