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
# preprocessing = false
preprocessing = true
# preprocess2_limit = 1
preprocess2_limit = 50
preprocess2_check_limit = 100_000
preprocess2_equal_poles = false
preprocess2_equal_poles = true
# preprocess3_limit = 1
preprocess3_limit = 10
preprocess3_check_limit = 100_000
preprocess3_equal_meters = false
# preprocess3_equal_meters = true
preprocess_repeats = 1
# preprocess_repeats = 3
# benchmarkTime = false
benchmarkTime = true

g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats,
    scoring_function=scoring_function,
    preprocessing=preprocessing,
    preprocess2_limit=preprocess2_limit,
    preprocess2_check_limit=preprocess2_check_limit,
    preprocess2_equal_poles=preprocess2_equal_poles,
    preprocess3_limit=preprocess3_limit,
    preprocess3_check_limit=preprocess3_check_limit,
    preprocess3_equal_meters=preprocess3_equal_meters,
    preprocess_repeats=preprocess_repeats)

#region solveSetCoveringProblem
if benchmarkTime
    sim_time = @elapsed begin
        Profile.clear()
        Profile.init()
        @profile begin
            local g_local = deepcopy(g)
            sCH.solveSetCoveringProblem!(g_local)
        end
    end
    global graph = g_local  # store the result (if needed)
else
    sim_time = "not_benchmarked"
end
#endregion solveSetCoveringProblem

@unpack poles_used, meters_covered, m = graph
@test meters_covered == m

myprintln(true, "************************")
println("Simulation Time[s]: $(sim_time)")
myprintln(true, "value=$(poles_used)")
myprintln(true, "scoring_function=$(scoring_function)")
myprintln(true, "testCase=$(testCase)")
myprintln(true, "cleanupRepeats=$(cleanupRepeats)")
myprintln(true, "************************")

@unpack poles_used, meters_covered, m, p, Pprime, Premaining, Pdiscarded, A_m2p, Aadj_m2p, A_m2p_remaining, Aadj_m2p_remaining, A_p2m, Aadj_p2m = graph;

@unpack preprocess1_steps, preprocess2_steps, preprocess3_steps = graph;
HF.myprintln(true, "preprocess1_steps=$(preprocess1_steps)")
HF.myprintln(true, "preprocess2_steps=$(preprocess2_steps)")
HF.myprintln(true, "preprocess3_steps=$(preprocess3_steps)")

@test meters_covered == m
open("profile_summary.txt", "w") do io
    Profile.print(io; format=:flat, sortedby=:count)
end