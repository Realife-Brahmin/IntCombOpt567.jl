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

using Combinatorics
using Parameters
using SparseArrays

include("./helperFunctions.jl")
import .helperFunctions as HF

function initializeGraph(filepath::String;
    maxiter::Int = 100000,
    cleanupRepeats::Int = 10,
    scoring_function::String = "greedy",
    preprocessing=false)
    # Initialize arrays to store row and column indices for A_m2p_remaining and A_p2m_uncovered
    rows_A = Int[]
    cols_A = Int[]
    rows_AT = Int[]
    cols_AT = Int[]

    # Initialize dictionaries to count degrees
    deg_m_uncovered = Dict{Int, Int}()  # Degree of each meter not currently in Mprime
    deg_p_remaining = Dict{Int, Int}()  # Degree of each pole (column j)
    deg_m2p_remaining = Dict{Int, Int}()   # Degree of each meter (row i) wrt poles NOT in Pprime
    deg_m2p = Dict{Int, Int}()  # Degree of each meter (row i) wrt poles part of Pprime

    # Read the file and parse the rows
    for line in readlines(filepath)
        i, j = parse.(Int, split(line))
        push!(rows_A, i)
        push!(cols_A, j)
        push!(rows_AT, j)  # Reverse for A_p2m_uncovered
        push!(cols_AT, i)

        # Update degree counts
        deg_m2p_remaining[i] = get(deg_m2p_remaining, i, 0) + 1
        deg_m_uncovered[i] = deg_m2p_remaining[i]
        deg_m2p[i] = 0
        deg_p_remaining[j] = get(deg_p_remaining, j, 0) + 1
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
    Aadj_m2p_remaining = Dict{Int,Vector{Int}}()
    Aadj_p2m_ref = Dict{Int,Vector{Int}}()
    Aadj_p2m = Dict{Int,Vector{Int}}()
    Aadj_p2m_uncovered = Dict{Int,Vector{Int}}()

    # Initialize adjacency lists
    for i in 1:max_i
        Aadj_m2p[i] = Vector{Int}()
        Aadj_m2p_ref[i] = Vector{Int}()
        Aadj_m2p_remaining[i] = Vector{Int}()
    end

    for j in 1:max_j
        Aadj_p2m_ref[j] = Vector{Int}()
        Aadj_p2m[j] = Vector{Int}()
        Aadj_p2m_uncovered[j] = Vector{Int}()
    end

    # Populate adjacency lists
    for (i, j) in zip(rows_A, cols_A)
        push!(Aadj_m2p_ref[i], j)  # Add pole j to the list of poles for meter i
        push!(Aadj_m2p_remaining[i], j)  # Add pole j to the list of remaining poles for meter i
        push!(Aadj_p2m_ref[j], i)  # Add meter i to the list of meters for pole j
        push!(Aadj_p2m_uncovered[j], i)  # Add meter i to the list of uncovered meters for pole j
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
        :Aadj_m2p_remaining => Aadj_m2p_remaining,

        :A_p2m => A_p2m,
        :Aadj_p2m => Aadj_p2m,
        :Aadj_p2m_ref => Aadj_p2m_ref,
        :A_p2m_uncovered => A_p2m_uncovered,
        :Aadj_p2m_uncovered => Aadj_p2m_uncovered,

        :cleanupDoneLastIter => false,
        :cleanupRepeats => cleanupRepeats,
        :cleanupUsefulLastIter => false,

        :preprocessing => preprocessing,

        :deg_p_remaining => deg_p_remaining,
        :deg_m_uncovered => deg_m_uncovered,
        :deg_m2p => deg_m2p,
        :deg_m2p_remaining => deg_m2p_remaining,

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

        :preprocess1_steps => 0,
        :preprocess2_steps => 0,
        :preprocess3_steps => 0,
    )

    return graphState
end

function chooseNextPole(graphState; verbose::Bool = false)

    @unpack meters_covered, m = graphState
    if meters_covered == m # Takes care of the case where preprocessing already covered all meters
        HF.myprintln(true, "All meters already covered. No pole selected.")
        return -1
    end

    @unpack deg_p_remaining, scoring_function = graphState
    if isempty(deg_p_remaining)
        error("No unused poles available.")
    end

    if scoring_function == "score1" || scoring_function == "score2"
        @unpack deg_m_uncovered, Aadj_m2p_remaining = graphState
        k = parse(Int, last(scoring_function))  # Extract the number from "score1" or "score2"
        scoreDict = compute_score(deg_m_uncovered, deg_p_remaining, Aadj_m2p_remaining, k=k, verbose=verbose)
    elseif scoring_function == "greedy"
        scoreDict = deg_p_remaining 
    else
        @error("Invalid scoring function: $scoring_function")
    end

    j_candidate = HF.argmax_smallestkey(scoreDict)[1]  # Pole with the maximum degree
    # j_candidate = argmax(scoreDict)[1]  # Pole with the maximum degree

    return j_candidate
end

function addPole!(graphState, j;
    verbose = false)

    if j == -1 # Takes care of the case where preprocessing already covered all meters
        HF.myprintln(true, "No pole selected for addition! Returning graph state unchanged.")
        return graphState
    end

    HF.myprintln(verbose, "Pole $j to be added")

    @unpack Aadj_p2m_ref = graphState
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")

    @unpack Pprime, Premaining, poles_used = graphState
    union!(Pprime, j)  # Add pole j to Pprime
    setdiff!(Premaining, j)  # Selected pole no longer on market
    poles_used = length(Pprime)
    @pack! graphState = Premaining, Pprime, poles_used # Update the graph state

    @unpack deg_p_remaining = graphState
    delete!(deg_p_remaining, j) # Selected pole no longer on market
    @pack! graphState = deg_p_remaining

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j

        @unpack A_m2p, Aadj_m2p, deg_m2p, A_m2p_remaining, deg_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered, Aadj_p2m_uncovered, deg_m_uncovered = graphState

        delete!(deg_m_uncovered, i)  # meter i now no longer uncovered (if previously uncovered at all)

        A_m2p[i, j] = 1
        push!(Aadj_m2p[i], j)  # Add pole j to the adjacency list of meter i
        A_m2p_remaining[i, j] = 0 # Remove pole j from the remaining poles for meter i
        setdiff!(Aadj_m2p_remaining[i], j)  # Remove pole j from the adjacency list of meter i
        deg_m2p_remaining[i] -= 1 # one less remaining pole which can cover meter i
        deg_m2p[i] += 1 # another pole which now covers meter i

        A_p2m[j, i] = 1
        push!(Aadj_p2m[j], i)  # Add meter i to the adjacency list of pole j
        A_p2m_uncovered[j, i] = 0 # Remove meter i from the uncovered meters for pole j

        # @show Aadj_p2m_uncovered
        setdiff!(Aadj_p2m_uncovered[j], i)  # Remove meter i from the adjacency list of pole j
        # @show Aadj_p2m_uncovered

        @pack! graphState = A_m2p, Aadj_m2p, deg_m2p, A_m2p_remaining, deg_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered, Aadj_p2m_uncovered, deg_m_uncovered

        if !(i ∈ graphState[:Mprime]) # i.e. a previously uncovered meter is being covered by pole j
            # delete!(deg_m_uncovered, i)  # meter i now no longer uncovered

            # All poles covering meter i still available in the market (excluding newly added pole j, already removed from Premaining)
            @unpack Aadj_m2p_ref = graphState
            other_remaining_poles_also_covering_i = intersect(Aadj_m2p_ref[i], Premaining)

            for j_other in other_remaining_poles_also_covering_i
                @unpack A_p2m_uncovered, Aadj_p2m_uncovered, deg_p_remaining = graphState

                HF.myprintln(verbose, "Pole $j_other also covers meter $i, which is now covered by pole $j")
                # @show Aadj_p2m_uncovered
                HF.myprintln(verbose, "Pole $j_other's previous list of uncovered meters was: $(Aadj_p2m_uncovered[j_other])")
                HF.myprintln(verbose, "Pole $j_other's previous degree was: $(deg_p_remaining[j_other])")
                A_p2m_uncovered[j_other, i] = 0  # Remove meter i from the uncovered meters for pole j_other
                setdiff!(Aadj_p2m_uncovered[j_other], i)  # Remove meter i from the adjacency list of pole j_other
                deg_p_remaining[j_other] -= 1  # Now that a meter is covered without pole j_other's help, its degree is deducted by 1

                HF.myprintln(verbose, "Pole $j_other's new list of uncovered meters is: $(Aadj_p2m_uncovered[j_other])")
                HF.myprintln(verbose, "Pole $j_other's new degree is: $(deg_p_remaining[j_other])")
                @pack! graphState = A_p2m_uncovered, Aadj_p2m_uncovered, deg_p_remaining
            end

        end
    end

    @unpack Mprime, meters_covered = graphState
    union!(Mprime, meters_covered_by_j)  # Add these meters to Mprime
    # Modifying Mprime only now as we need to check if the meters are already in Mprime

    meters_covered = length(Mprime) 
    # Update the graph state
    @pack! graphState = Mprime, meters_covered
    return graphState
end

function compute_score(deg_m_uncovered, deg_p_remaining, Aadj_m2p_remaining; verbose::Bool = false, k=1)

    htc_threshold = minimum(values(deg_m_uncovered))  # Hard-to-cover threshold (minimum degree of uncovered meters)
    HF.myprintln(verbose, "Hard-to-cover threshold: $htc_threshold")
    scoreDict = Dict{Int, Int}()
    
    if k == 1
        for i ∈ keys(deg_m_uncovered)
            if deg_m_uncovered[i] == htc_threshold
                # HF.myprintln(true, "Meter $i is a hard to cover meter")
                for j ∈ Aadj_m2p_remaining[i]
                    # HF.myprintln(true, "Meter $i is covered by pole $j")
                    scoreDict[j] = deg_p_remaining[j]
                end
            end
        end
    elseif k == 2
        for i ∈ keys(deg_m_uncovered)
            if deg_m_uncovered[i] == htc_threshold
                # HF.myprintln(true, "Meter $i is a hard to cover meter")
                for j ∈ Aadj_m2p_remaining[i]
                    # HF.myprintln(true, "Meter $i is covered by pole $j")
                    scoreDict[j] = get(scoreDict, j, 0) + deg_p_remaining[j]
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
    @unpack Premaining, deg_p_remaining, deg_m2p, deg_m2p_remaining, Mprime, Pprime  = graphState

    if j ∉ Pprime
        error("Attempting to remove a pole that is not in P′.")
        return
    end

    setdiff!(Pprime, j)  # Remove pole j from Pprime
    # But j will not be placed back in Premaining (by design)

    HF.myprintln(verbose, "Pole $j to be removed")

    @unpack Aadj_p2m_ref = graphState
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")
    
    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j

        @unpack A_m2p, Aadj_m2p, A_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered = graphState;
        A_m2p[i, j] = 0  # Remove pole j from the sparse matrix A_m2p
        setdiff!(Aadj_m2p[i], j)  # Remove pole j from the adjacency list of meter i
        A_m2p_remaining[i, j] = 1  # Add pole j back to the sparse matrix A_m2p_remaining
        push!(Aadj_m2p_remaining[i], j)  # Add pole j back to the adjacency list of meter i
        A_p2m[j, i] = 0  # Remove meter i from the sparse matrix A_p2m
        setdiff!(Aadj_p2m[j], i)  # Remove meter i from the adjacency list of pole j
        A_p2m_uncovered[j, i] = 1  # Add meter i back to the sparse matrix A_p2m_uncovered
        @pack! graphState = A_m2p, Aadj_m2p, A_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered

        deg_m2p_remaining[i] += 1  # Update degrees for meters covered by pole j
        deg_m2p[i] -= 1
        if deg_m2p[i] == 0
            @error("Removal of pole $j leaves meter $i uncovered! Why was this selected?")
        end
    end

    poles_used = length(Pprime)
    meters_covered = length(Mprime)
    # Update the graph state
    @pack! graphState = Mprime, Pprime, poles_used, deg_p_remaining, deg_m2p, deg_m2p_remaining, meters_covered
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

        if graphState[:preprocessing]
            preprocess1!(graphState; verbose=verbose)  # Preprocess the graph to find singleton meters
            preprocess2!(graphState; verbose=true)  # Preprocess the graph to find dominating poles
        end
        j = chooseNextPole(graphState)  # Choose the next pole
        addPole!(graphState, j, verbose=verbose)  # Select the pole and update the graph state

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

    sanitize_data_structures!(graphState; verbose=verbose)  # Clean up the matrices for the final output

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
    @unpack deg_m2p, Aadj_p2m_ref = graphState
    meters_covered_by_j = Aadj_p2m_ref[j]  # Find all meters covered by pole j
    
    for i in meters_covered_by_j
        if deg_m2p[i] == 1  # If meter i is only covered by pole j
            HF.myprintln(verbose, "Pole $j is NOT redundant as it covers meter $i exclusively")
            return false
        end
    end

    return true
end

function preprocess1!(graphState;
    verbose::Bool = false)

    keepPP1Running = true
    while keepPP1Running
        @unpack deg_m_uncovered, Aadj_m2p_remaining = graphState
        graph_mutated = false
        for i ∈ keys(deg_m_uncovered)
            if deg_m_uncovered[i] == 1  # If meter i is only covered by one pole
                HF.myprintln(verbose, "Meter $i is a singleton meter")
                if length(Aadj_m2p_remaining[i]) != 1
                    @error("Meter $i is a singleton meter but has more than one pole covering it!")
                end
                j = Aadj_m2p_remaining[i][1]  # The only pole covering meter i
                HF.myprintln(verbose, "Pole $j covers meter $i exclusively")
                addPole!(graphState, j; verbose=verbose)  # Add the pole to Pprime
                graphState[:preprocess1_steps] += 1  # Increment the number of steps taken in preprocess1
                graph_mutated = true
            end

            if graph_mutated # Adding a pole mutates the graph, mutating its fields like deg_m_uncovered, etc. in the process. We can no longer iterate over the already unpacked (unmutated) fields.
                break;
            end
        end

        if !graph_mutated
            HF.myprintln(verbose, "No more singleton meters found")
            keepPP1Running = false  # No more meters to process
        end
    end

    return graphState
end

function preprocess2!(graphState;
    verbose::Bool=false)

    keepPP2Running = true
    while keepPP2Running
        @unpack deg_p_remaining, Aadj_p2m_uncovered = graphState
        graph_mutated = false
        for (j1, j2) in combinations(collect(keys(deg_p_remaining)), 2)
            # HF.myprintln(verbose, "Comparing poles $j1 and $j2")
            # HF.myprintln(verbose, "Pole $j1 has degree $(deg_p_remaining[j1])")
            # HF.myprintln(verbose, "Pole $j2 has degree $(deg_p_remaining[j2])")
            # HF.myprintln(verbose, "Pole $j1 covers meters: $(Aadj_p2m_uncovered[j1]), degree = $(deg_p_remaining[j1])")
            # HF.myprintln(verbose, "Pole $j2 covers meters: $(Aadj_p2m_uncovered[j2]), degree = $(deg_p_remaining[j2])")
            # if deg_p_remaining[j1] != deg_p_remaining[j2]
            if deg_p_remaining[j1] < deg_p_remaining[j2]
                j_small, j_big = j1, j2
            else
                j_small, j_big = j2, j1
            end

            # Check dominance: j_big covers all meters covered by j_small
            if issubset(Aadj_p2m_uncovered[j_small], Aadj_p2m_uncovered[j_big])
                HF.myprintln(verbose, "Pole $(j_big) is dominant over pole $(j_small).")
                
                discardPole!(graphState, j_small; verbose=verbose)  # Remove the smaller pole from the graph
                graph_mutated = true
                graphState[:preprocess2_steps] += 1
            else
                # HF.myprintln(verbose, "Pole $(j_big) is NOT dominant over pole $(j_small).")
            end
            # end

            if graph_mutated # Discarding a pole mutates the graph, mutating its fields like deg_p_remaining, etc. in the process. We can no longer iterate over the already unpacked (unmutated) fields.
                break
            end
        end

        if !graph_mutated
            HF.myprintln(verbose, "No more dominating poles found")
            keepPP2Running = false  # No more meters to process
        end
    end

    return graphState
end

function discardPole!(graphState, j;
    verbose::Bool = false) # discardPole! removes a 'remaining' pole from graph, unlike removePole! which removes a pole from Pprime. Currently meant for use in preprocess2!

    @unpack Premaining, deg_p_remaining = graphState

    if j ∉ Premaining
        error("Attempting to discard a pole that is not in Premaining.")
        return
    end

    HF.myprintln(verbose, "Pole $j to be discarded")

    setdiff!(Premaining, j)  # Remove pole j from Premaining
    delete!(deg_p_remaining, j) # Selected pole no longer on market

    @pack! graphState = Premaining, deg_p_remaining
    
    @unpack A_m2p_remaining, Aadj_m2p_remaining, deg_m_uncovered, deg_m2p_remaining, Aadj_p2m_uncovered = graphState

    meters_covered_by_j = Aadj_p2m_uncovered[j]  # Find all meters covered by pole j
    # Update degrees for meters covered by pole j
    # @show meters_covered_by_j
    for i in meters_covered_by_j

        A_m2p_remaining[i, j] = 0  # Disqualifying pole j as a 'remaining' pole for meter i
        setdiff!(Aadj_m2p_remaining[i], j)  # Disqualifying pole j as a 'remaining' pole for meter i
        
        deg_m2p_remaining[i] -= 1  # Pole j is no longer a remaining pole for meter i

        if !(i ∈ graphState[:Mprime])
            deg_m_uncovered[i] -= 1  # Pole j is no longer a remaining pole for meter i (if previously uncovered at all)
        end
    end

    # Update the graph state
    @pack! graphState = A_m2p_remaining, Aadj_m2p_remaining, deg_m_uncovered, deg_m2p_remaining

    return graphState
end

function sanitize_data_structures!(graphState;
    verbose::Bool = false)

    @unpack A_m2p, A_m2p_remaining, A_p2m, A_p2m_uncovered = graphState

    dropzeros!(A_m2p)  # Remove zero rows from A_m2p
    dropzeros!(A_m2p_remaining)  # Remove zero rows from A_m2p_remaining
    dropzeros!(A_p2m)  # Remove zero rows from A_p2m
    dropzeros!(A_p2m_uncovered)  # Remove zero rows from A_p2m_uncovered

    return graphState
end

end # module setCoveringHeuristics