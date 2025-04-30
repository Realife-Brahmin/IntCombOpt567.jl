# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

# testCase = "p4m6"
# testCase = "p6m6"
# testCase = "phase1"
testCase = "cap360"
cleanupRepeats = 1
# cleanupRepeats = 10
# cleanupRepeats = 30
# cleanupRepeats = 100
scoring_function = "greedy"
# scoring_function = "score1"
# scoring_function = "score2"
benchmarkTime = false
# benchmarkTime = true

g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats,
    scoring_function=scoring_function)

myprintln(true, "************************")

#region solveSetCoveringProblem
if benchmarkTime
    # just time the solve call
    sim_time = @belapsed begin
        g_local = deepcopy(g)
        sCH.solveSetCoveringProblem!(g_local)
    end
    println("Simulation Time: $(sim_time) s")
else
    sim_time = "not_benchmarked"
    println("Simulation Time: $(sim_time)")
end

# profile separately (only once)
@profile begin
    sCH.solveSetCoveringProblem!(g)
    global graph = g
end
#endregion solveSetCoveringProblem

@unpack poles_used, meters_covered, m = graph
@test meters_covered == m

myprintln(true, "value=$(poles_used)")
myprintln(true, "scoring_function=$(scoring_function)")
myprintln(true, "testCase=$(testCase)")
myprintln(true, "cleanupRepeats=$(cleanupRepeats)")
myprintln(true, "************************")

@unpack poles_used, meters_covered, m, A_m2p, A_m2p_remaining, A_p2m = graph;

poles_used_as_per_Am2p = length(unique(findnz(A_m2p)[2]))  # Extract column indices from non-zero entries
poles_not_used_as_per_Am2p_remaining = length(unique(findnz(A_m2p_remaining)[2]))  # Extract column indices from non-zero entries
poles_used_as_per_Am2p_remaining = p - poles_not_used_as_per_Am2p_remaining

HF.myprintln(true, "Poles used as per A_m2p: $(poles_used_as_per_Am2p)")
HF.myprintln(true, "Poles used as per A_m2p_remaining: $(poles_used_as_per_Am2p_remaining)")
HF.myprintln(true, "Poles used as per poles_used: $(poles_used)")


# open("profile_summary.txt", "w") do io
#     Profile.print(io; format=:flat, sortedby=:count)
# end