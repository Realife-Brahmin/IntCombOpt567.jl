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
    scoring_function::String = "greedy",
    preprocessing=false,
    preprocess2_limit=1,
    preprocess2_check_limit=10,
    preprocess2_equal_poles=false,
    preprocess3_limit=1,
    preprocess3_check_limit=10,
    preprocess3_equal_meters=false,
    preprocess_repeats=1
    )
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
    harder_to_cover_meter = Dict{Int, Int}()  # Should meter i be ignored, will map to the meter, covering which would also cover it

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

        harder_to_cover_meter[i] = -1
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

        :harder_to_cover_meter => harder_to_cover_meter,

        :k => 0,
        :M0 => M,
        :M => M,
        :m0 => m,
        :m => m,
        :maxiter => maxiter,
        :Mprime => Mprime,
        :Mignored => Set{Int}(),  # Set of ignored meters (initially empty)
        :P => P,
        :Premaining => Premaining,
        :Pdiscarded => Set{Int}(),  # Set of discarded meters (initially empty)
        :p => p,
        :Pprime => Pprime,
        :poles_used => 0,
        :meters_covered => 0,
        :scoring_function => scoring_function,

        :preprocess1_steps => 0,
        :preprocess2_steps => 0,
        :preprocess2_limit => preprocess2_limit,
        :preprocess2_check_limit => preprocess2_check_limit,
        :preprocess2_equal_poles => preprocess2_equal_poles,
        :preprocess3_steps => 0,
        :preprocess3_limit => preprocess3_limit,
        :preprocess3_check_limit => preprocess3_check_limit,
        :preprocess3_equal_meters => preprocess3_equal_meters,
        :preprocess_repeats => preprocess_repeats,
    )

    return graphState
end

function chooseNextPole(graphState; verbose::Bool = false)

    @unpack m = graphState
    # @show graphState[:meters_covered]    
    if graphState[:meters_covered] == m # Takes care of the case where preprocessing already covered all meters
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

    # j_candidate = HF.argmax_smallestkey(scoreDict)[1]  # Pole with the maximum degree
    j_candidate = HF.argmax(scoreDict)[1]  # Pole with the maximum degree

    # j_candidate = argmax(scoreDict)[1]  # Pole with the maximum degree

    return j_candidate
end

function addPole!(graphState, j;
    verbose=false)

    if j == -1 # Takes care of the case where preprocessing already covered all meters
        HF.myprintln(true, "No pole selected for addition! Returning graph state unchanged.")
        return graphState
    end

    HF.myprintln(verbose, "Pole $j to be added")

    @unpack Aadj_p2m_ref = graphState
    meters_covered_by_j = deepcopy(Aadj_p2m_ref[j])  # Find all meters covered by pole j
    meters_covered_by_j = setdiff(meters_covered_by_j, graphState[:Mignored])
    # @show typeof(meters_covered_by_j)
    # @assert pointer_from_objref(meters_covered_by_j) != pointer_from_objref(Aadj_p2m_uncovered[j])

    @unpack Pprime, Premaining, poles_used = graphState
    union!(Pprime, j)  # Add pole j to Pprime
    setdiff!(Premaining, j)  # Selected pole no longer on market
    # @show j ∈ Premaining
    poles_used = length(Pprime)
    @pack! graphState = Premaining, Pprime, poles_used # Update the graph state

    @unpack deg_p_remaining = graphState
    delete!(deg_p_remaining, j) # Selected pole no longer on market
    # @show j ∈ keys(deg_p_remaining)
    @pack! graphState = deg_p_remaining

    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j

        @unpack A_m2p, Aadj_m2p, deg_m2p, A_m2p_remaining, deg_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered, Aadj_p2m_uncovered, deg_m_uncovered = graphState

        delete!(deg_m_uncovered, i)  # meter i now no longer uncovered (if previously uncovered at all)

        A_m2p[i, j] = 1
        push!(Aadj_m2p[i], j)  # Add pole j to the adjacency list of meter i
        A_m2p_remaining[i, j] = 0 # Remove pole j from the remaining poles for meter i
        setdiff!(Aadj_m2p_remaining[i], j)  # Remove pole j from the adjacency list of meter i
        deg_m2p[i] += 1 # another pole which now covers meter i
        deg_m2p_remaining[i] -= 1 # one less remaining pole which can cover meter i

        A_p2m[j, i] = 1
        push!(Aadj_p2m[j], i)  # Add meter i to the adjacency list of pole j
        A_p2m_uncovered[j, i] = 0 # Remove meter i from the uncovered meters for pole j

        # @show Aadj_p2m_uncovered
        setdiff!(Aadj_p2m_uncovered[j], i)  # Remove meter i from the adjacency list of pole j
        # @show Aadj_p2m_uncovered

        @pack! graphState = A_m2p, Aadj_m2p, deg_m2p, A_m2p_remaining, deg_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered, Aadj_p2m_uncovered, deg_m_uncovered

        if !(i ∈ graphState[:Mprime]) # i.e. a previously uncovered meter is being covered by pole j
            @unpack Aadj_m2p_remaining = graphState
            other_remaining_poles_also_covering_i = deepcopy(Aadj_m2p_remaining[i])  # Find all remaining poles covering meter i
            for j_other in other_remaining_poles_also_covering_i
                @unpack A_p2m_uncovered, Aadj_p2m_uncovered, deg_p_remaining = graphState

                # HF.myprintln(verbose, "Pole $j_other also covers meter $i, which is now covered by pole $j")
                # @show Aadj_p2m_uncovered
                # HF.myprintln(verbose, "Pole $j_other's previous list of uncovered meters was: $(Aadj_p2m_uncovered[j_other])")
                # HF.myprintln(verbose, "Pole $j_other's previous degree was: $(deg_p_remaining[j_other])")
                A_p2m_uncovered[j_other, i] = 0  # Remove meter i from the uncovered meters for pole j_other
                setdiff!(Aadj_p2m_uncovered[j_other], i)  # Remove meter i from the adjacency list of pole j_other
                deg_p_remaining[j_other] -= 1  # Now that a meter is covered without pole j_other's help, its degree is deducted by 1

                # HF.myprintln(verbose, "Pole $j_other's new list of uncovered meters is: $(Aadj_p2m_uncovered[j_other])")
                # HF.myprintln(verbose, "Pole $j_other's new degree is: $(deg_p_remaining[j_other])")
                @pack! graphState = A_p2m_uncovered, Aadj_p2m_uncovered, deg_p_remaining
            end

        end
    end

    @unpack Mprime = graphState
    # @show Mprime
    # @show meters_covered_by_j
    union!(Mprime, meters_covered_by_j)  # Add these meters to Mprime
    # Modifying Mprime only now as we need to check if the meters are already in Mprime

    graphState[:meters_covered] = length(Mprime) 
    # Update the graph state
    @pack! graphState = Mprime
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

    HF.myprintln(verbose, "Removing and discarding pole $j")
    union!(graphState[:Pdiscarded], j)  # Add pole j to Pdiscarded
    
    # HF.myprintln(verbose, "Pole $j to be removed")

    @unpack Aadj_p2m_ref = graphState
    meters_covered_by_j = deepcopy(Aadj_p2m_ref[j])  # Find all meters covered by pole j. It is okay to use the reference matrix here, as all meters in it must be in Mprime.
    meters_covered_by_j = setdiff(meters_covered_by_j, graphState[:Mignored])
    HF.myprintln(verbose, "Pole $j covers meters:  $(meters_covered_by_j)")
    
    # Update degrees for meters covered by pole j
    for i in meters_covered_by_j

        @unpack A_m2p, Aadj_m2p, A_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered = graphState;
        A_m2p[i, j] = 0  # Remove pole j from the sparse matrix A_m2p
        setdiff!(Aadj_m2p[i], j)  # Remove pole j from the adjacency list of meter i
        # A_m2p_remaining[i, j] = 1  # Add pole j back to the sparse matrix A_m2p_remaining
        # push!(Aadj_m2p_remaining[i], j)  # Add pole j back to the adjacency list of meter i
        A_p2m[j, i] = 0  # Remove meter i from the sparse matrix A_p2m
        setdiff!(Aadj_p2m[j], i)  # Remove meter i from the adjacency list of pole j
        # A_p2m_uncovered[j, i] = 1  # Add meter i back to the sparse matrix A_p2m_uncovered
        @pack! graphState = A_m2p, Aadj_m2p, A_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m, A_p2m_uncovered

        deg_m2p_remaining[i] += 1  # Update degrees for meters covered by pole j
        deg_m2p[i] -= 1
        if deg_m2p[i] == 0
            @error("Removal of pole $j leaves meter $i uncovered! Why was this selected?")
        end
    end

    poles_used = length(Pprime)
    graphState[:meters_covered] = length(Mprime)
    # Update the graph state
    @pack! graphState = Mprime, Pprime, poles_used, deg_p_remaining, deg_m2p, deg_m2p_remaining
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
            preprocess_steps_this_iter = 0
            while preprocess_steps_this_iter < graphState[:preprocess_repeats]
                preprocess1!(graphState; verbose=true)  # Preprocess the graph to find singleton meters
                preprocess2!(graphState; verbose=true)  # Preprocess the graph to find dominating poles
                preprocess3!(graphState; verbose=true)  # Preprocess the graph to find hard-to-cover meters
                preprocess_steps_this_iter += 1
            end
        end
        j = chooseNextPole(graphState,)  # Choose the next pole
        addPole!(graphState, j, verbose=true)  # Select the pole and update the graph state
        # @show graphState[:meters_covered]
        @pack! graphState = k # k-th iteration completed, so saving it

        cleanupAttempt = false
        if graphState[:meters_covered] == m
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
    meters_covered_by_j = setdiff(meters_covered_by_j, graphState[:Mignored])
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
                # HF.myprintln(verbose, "Meter $i is a singleton meter")
                if length(Aadj_m2p_remaining[i]) != 1
                    @error("Meter $i is a singleton meter but has more than one pole covering it!")
                end
                j = Aadj_m2p_remaining[i][1]  # The only pole covering meter i
                # HF.myprintln(verbose, "Pole $j covers meter $i exclusively, so adding it to P′")
                addPole!(graphState, j; verbose=verbose)  # Add the pole to Pprime
                # @show graphState[:meters_covered]
                graphState[:preprocess1_steps] += 1  # Increment the number of steps taken in preprocess1
                graph_mutated = true
            end

            if graph_mutated # Adding a pole mutates the graph, mutating its fields like deg_m_uncovered, etc. in the process. We can no longer iterate over the already unpacked (unmutated) fields.
                break;
            end
        end

        if !graph_mutated
            # HF.myprintln(verbose, "No more singleton meters found")
            keepPP1Running = false  # No more meters to process
        end
    end

    return graphState
end

function preprocess2!(graphState; verbose::Bool=false)
    @unpack preprocess2_limit, preprocess2_check_limit, preprocess2_equal_poles = graphState
    preprocess2_steps_this_iter = 0
    preprocess2_check_steps_this_iter = 0
    keepPP2Running = true

    while keepPP2Running
        @unpack Premaining, deg_p_remaining, Aadj_p2m_uncovered = graphState
        # HF.myprintln(true, "PP2: Current set of remaining poles: $(keys(deg_p_remaining))")

        poles = collect(Premaining)
        poles_to_remove = Int[]
        graph_mutated = false

        # Sweep through all pairs, skipping any already marked for deletion
        # for (j1, j2) in combinations(poles, 2)
        # for idx1 ∈ eachindex(poles), idx2 = idx1+1:length(poles)
        #     j1, j2 = poles[idx1], poles[idx2]
        #     if j1 ∈ poles_to_remove || j2 ∈ poles_to_remove
        #         continue
        #     end

        #     if !preprocess2_equal_poles && deg_p_remaining[j1] == deg_p_remaining[j2]
        #         continue  # Skip if the poles are equal and we want to ignore them
        #     end
        #     # Determine smaller and larger degree poles

        #     if deg_p_remaining[j1] == 0
        #         HF.myprintln(verbose, "Pole $j1 has no remaining degree. Forwarding for removal")
        #         push!(poles_to_remove, j1)
        #         preprocess2_steps_this_iter += 1
        #         continue
        #     end
        #     if deg_p_remaining[j2] == 0
        #         HF.myprintln(verbose, "Pole $j2 has no remaining degree. Forwarding for removal")
        #         push!(poles_to_remove, j2)
        #         preprocess2_steps_this_iter += 1
        #         continue
        #     end

        #     if deg_p_remaining[j1] < deg_p_remaining[j2]
        #         j_small, j_big = j1, j2
        #     else
        #         j_small, j_big = j2, j1
        #     end

        #     meters_small, meters_big = Set(Aadj_p2m_uncovered[j_small]), Set(Aadj_p2m_uncovered[j_big])
        #     # Check if j_big dominates j_small
        #     if issubset(meters_small, meters_big)
        #         HF.myprintln(verbose, "Pole $j_big dominates pole $j_small. Forwarding pole $j_small for removal")
        #         push!(poles_to_remove, j_small)
        #         preprocess2_steps_this_iter += 1
        #     end

        #     preprocess2_check_steps_this_iter += 1
        #     if preprocess2_check_steps_this_iter >= preprocess2_check_limit
        #         break  # Stop if the limit is reached
        #     end
        #     if preprocess2_steps_this_iter >= preprocess2_limit
        #         break  # Stop if the limit is reached
        #     end
        # end

        for idx1 in eachindex(poles)
            j1 = poles[idx1]
            if j1 in poles_to_remove
                continue  # Skip this j1 entirely
            end

            # Check and remove j1 or j2 if degree is 0
            if deg_p_remaining[j1] == 0
                HF.myprintln(verbose, "Pole $j1 has no remaining degree. Forwarding for removal")
                push!(poles_to_remove, j1)
                preprocess2_steps_this_iter += 1
                # break  # j1 is no longer valid for other j2s
                continue
            end
            
            for idx2 in (idx1+1):length(poles)
                j2 = poles[idx2]
                if j2 in poles_to_remove
                    continue  # Skip this j2 only
                end

                if !preprocess2_equal_poles && deg_p_remaining[j1] == deg_p_remaining[j2]
                    continue
                end

                if deg_p_remaining[j2] == 0
                    HF.myprintln(verbose, "Pole $j2 has no remaining degree. Forwarding for removal")
                    push!(poles_to_remove, j2)
                    preprocess2_steps_this_iter += 1
                    continue  # j2 removed, try next
                end

                # Determine j_small and j_big
                if deg_p_remaining[j2] <= deg_p_remaining[j2]
                    j_small, j_big = j1, j2
                else
                    j_small, j_big = j2, j1
                end

                meters_small = Set(Aadj_p2m_uncovered[j_small])
                meters_big = Set(Aadj_p2m_uncovered[j_big])

                if issubset(meters_small, meters_big)
                    # HF.myprintln(verbose, "Pole $j_big dominates pole $j_small. Forwarding pole $j_small for removal")
                    push!(poles_to_remove, j_small)
                    preprocess2_steps_this_iter += 1

                    if j_small == j1
                        break  # no point testing j1 further
                    else
                        continue  # skip this j2, continue for j1
                    end
                end

                preprocess2_check_steps_this_iter += 1
                if preprocess2_check_steps_this_iter >= preprocess2_check_limit || preprocess2_steps_this_iter >= preprocess2_limit
                    break
                end
            end

            # Optional: check global limits again here to fully break both loops
            if preprocess2_check_steps_this_iter >= preprocess2_check_limit ||preprocess2_steps_this_iter >= preprocess2_limit
                break
            end
        end

        # Remove dominated poles in one batch
        if !isempty(poles_to_remove)
            graph_mutated = true
            for j in poles_to_remove
                discardPole!(graphState, j; verbose=verbose)
                graphState[:preprocess2_steps] += 1
            end
        end

        keepPP2Running = false
    end

    return graphState
end

function ignoreMeter!(graphState, i::Int; verbose::Bool=false)
    @unpack A_m2p, Aadj_m2p, Aadj_m2p_ref, A_m2p_remaining, Aadj_m2p_remaining = graphState
    @unpack A_p2m_uncovered, Aadj_p2m_uncovered, A_p2m, Aadj_p2m = graphState
    @unpack deg_m2p, deg_m2p_remaining, deg_m_uncovered = graphState
    @unpack deg_p_remaining = graphState

    push!(graphState[:Mignored], i)

    # --- For every pole j connected to meter i ---
    poles_covering_i = deepcopy(Aadj_m2p_ref[i])
    poles_covering_i = setdiff(poles_covering_i, graphState[:Pdiscarded])  # Remove ignored poles from the list of covering poles
    for j in poles_covering_i
        # Sparse matrix cleanup
        A_m2p[i, j] = 0
        A_m2p_remaining[i, j] = 0
        A_p2m[j, i] = 0
        A_p2m_uncovered[j, i] = 0

        # Adjacency list cleanup
        setdiff!(Aadj_p2m[j], i)
        setdiff!(Aadj_p2m_uncovered[j], i)

        deg_p_remaining[j] -= 1  # Update pole's degree
        # # Update pole's remaining degree if meter i is still uncovered
        # if i ∉ graphState[:Mprime] # which should always be the case
        #     deg_p_remaining[j] = get(deg_p_remaining, j, 1) - 1
        # end
    end

    # --- Delete meter i’s entries ---
    delete!(Aadj_m2p, i)
    # delete!(Aadj_m2p_ref, i)
    delete!(Aadj_m2p_remaining, i)
    delete!(deg_m2p, i)
    delete!(deg_m2p_remaining, i)
    delete!(deg_m_uncovered, i)

    # --- Remove from M and Mprime ---
    setdiff!(graphState[:M], i)
    # --- Update m ---
    graphState[:m] = length(graphState[:M])

    HF.myprintln(verbose, "Ignoring meter $i")
    # HF.myprintln(verbose, "Remaining meters: $(graphState[:M])")

    return graphState
end

# function preprocess3!(graphState; verbose::Bool=false)
#     @unpack preprocess3_limit, preprocess3_check_limit, preprocess3_equal_meters = graphState
#     preprocess3_steps_this_iter = 0
#     preprocess3_check_steps_this_iter = 0

#     @unpack harder_to_cover_meter, deg_m_uncovered, Aadj_m2p_remaining, deg_m2p_remaining = graphState
#     meters = collect(keys(deg_m_uncovered))  # All uncovered meters
#     meters_to_ignore = Int[]

#     # Sweep through pairs of meters
#     # for (i1, i2) in combinations(meters, 2)
#     for idx1 ∈ eachindex(meters), idx2 = idx1+1:length(meters)
#         i1, i2 = meters[idx1], meters[idx2]
#         if i1 ∈ meters_to_ignore || i2 ∈ meters_to_ignore
#             continue
#         end

#         poles1 = Aadj_m2p_remaining[i1]
#         poles2 = Aadj_m2p_remaining[i2]

#         if !preprocess3_equal_meters && deg_m2p_remaining[i1] == deg_m2p_remaining[i2]
#             continue
#         end

#         # Decide which meter is easier to cover
#         if length(poles1) < length(poles2)
#             i_easy, i_hard = i2, i1
#             poles_more, poles_less = Set(poles2), Set(poles1)
#         else
#             i_easy, i_hard = i1, i2
#             poles_more, poles_less = Set(poles1), Set(poles2)
#         end

#         # Check if i_easy is dominated by i_hard (i.e., all poles that cover i_easy also cover i_hard)
#         if issubset(poles_less, poles_more)
#             push!(meters_to_ignore, i_easy)
#             harder_to_cover_meter[i_easy] = i_hard
#             preprocess3_steps_this_iter += 1
#         end

#         preprocess3_check_steps_this_iter += 1
#         if preprocess3_check_steps_this_iter >= preprocess3_check_limit
#             break
#         end
#         if preprocess3_steps_this_iter >= preprocess3_limit
#             break
#         end
#     end

#     # Perform ignoreMeter! for all meters found
#     for i in meters_to_ignore
#         ignoreMeter!(graphState, i; verbose=verbose)
#         graphState[:preprocess3_steps] += 1
#     end

#     return graphState
# end
function preprocess3!(graphState; verbose::Bool=false)
    @unpack preprocess3_limit, preprocess3_check_limit, preprocess3_equal_meters = graphState
    preprocess3_steps_this_iter = 0
    preprocess3_check_steps_this_iter = 0
    keepPP3Running = true

    while keepPP3Running
        @unpack harder_to_cover_meter, deg_m_uncovered, Aadj_m2p_remaining, deg_m2p_remaining = graphState
        meters = collect(keys(deg_m_uncovered))
        meters_to_ignore = Int[]
        graph_mutated = false

        for idx1 in eachindex(meters)
            i1 = meters[idx1]
            if i1 in meters_to_ignore
                continue
            end

            # Check and remove i1 if it has no adjacent poles left
            if deg_m2p_remaining[i1] == 0
                HF.myprintln(verbose, "Meter $i1 has no remaining poles. Forwarding for ignore")
                push!(meters_to_ignore, i1)
                preprocess3_steps_this_iter += 1
                continue
            end

            for idx2 in (idx1+1):length(meters)
                i2 = meters[idx2]
                if i2 in meters_to_ignore
                    continue
                end

                if !preprocess3_equal_meters && deg_m2p_remaining[i1] == deg_m2p_remaining[i2]
                    continue
                end

                if deg_m2p_remaining[i2] == 0
                    HF.myprintln(verbose, "Meter $i2 has no remaining poles. Forwarding for ignore")
                    push!(meters_to_ignore, i2)
                    preprocess3_steps_this_iter += 1
                    continue
                end

                if deg_m2p_remaining[i1] <= deg_m2p_remaining[i2]
                    i_easy, i_hard = i1, i2
                else
                    i_easy, i_hard = i2, i1
                end

                poles_easy = Set(Aadj_m2p_remaining[i_easy])
                poles_hard = Set(Aadj_m2p_remaining[i_hard])

                if issubset(poles_easy, poles_hard)
                    HF.myprintln(verbose, "Meter $i_hard dominates meter $i_easy. Forwarding $i_easy for ignore")
                    push!(meters_to_ignore, i_easy)
                    harder_to_cover_meter[i_easy] = i_hard
                    preprocess3_steps_this_iter += 1

                    if i_easy == i1
                        break  # Stop considering i1
                    else
                        continue  # Move to next i2
                    end
                end

                preprocess3_check_steps_this_iter += 1
                if preprocess3_check_steps_this_iter >= preprocess3_check_limit || preprocess3_steps_this_iter >= preprocess3_limit
                    break
                end
            end

            if preprocess3_check_steps_this_iter >= preprocess3_check_limit || preprocess3_steps_this_iter >= preprocess3_limit
                break
            end
        end

        if !isempty(meters_to_ignore)
            graph_mutated = true
            for i in meters_to_ignore
                ignoreMeter!(graphState, i; verbose=verbose)
                graphState[:preprocess3_steps] += 1
            end
        end

        keepPP3Running = false
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

    HF.myprintln(verbose, "Discarding pole $j")
    union!(graphState[:Pdiscarded], j)  # Add pole j to Pdiscarded

    # HF.myprintln(verbose, "Pole $j to be discarded")

    setdiff!(Premaining, j)  # Remove pole j from Premaining
    # @show deg_p_remaining
    delete!(deg_p_remaining, j) # Selected pole no longer on market
    # @show deg_p_remaining
    @pack! graphState = Premaining, deg_p_remaining
    
    @unpack A_m2p_remaining, Aadj_m2p_remaining, deg_m_uncovered, deg_m2p_remaining, Aadj_p2m_uncovered = graphState

    # HF.myprintln(true, "Discarded pole $j would have covered these uncovered meters:  $(Aadj_p2m_uncovered[j])")
    # meters_covered_by_j = deepcopy(Aadj_p2m_uncovered[j])  
    meters_covered_by_j = graphState[:Aadj_p2m_ref][j]  # Find all meters covered by pole j. Double remove them if necessary
    meters_covered_by_j = setdiff(meters_covered_by_j, graphState[:Mignored])
    # HF.myprintln(true, "Removing pole $j as a candidate for these meters: $(meters_covered_by_j)")
    # Find all meters covered by pole j
    
    
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