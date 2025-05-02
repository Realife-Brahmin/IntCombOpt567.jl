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
preprocessing = false
preprocessing = true
benchmarkTime = false
# benchmarkTime = true

g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats,
    scoring_function=scoring_function,
    preprocessing=preprocessing)

#region solveSetCoveringProblem
if benchmarkTime
    # just time the solve call
    sim_time = @belapsed begin
        g_local = deepcopy(g)
        sCH.solveSetCoveringProblem!(g_local)
    end
else
    sim_time = "not_benchmarked"
end

# profile separately (only once)
@profile begin
    sCH.solveSetCoveringProblem!(g)
    global graph = g
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

poles_used_as_per_Am2p = length(unique(findnz(A_m2p)[2]))  # Extract column indices from non-zero entries
poles_not_used_as_per_Am2p_remaining_and_Pdiscarded = length(unique(findnz(A_m2p_remaining)[2])) + length(Pdiscarded)  # Extract column indices from non-zero entries

poles_used_as_per_Am2p_remaining_and_Pdiscarded = p - poles_not_used_as_per_Am2p_remaining_and_Pdiscarded
poles_used_as_per_Aadj_m2p_remaining_and_Pdiscarded = p - length(unique(vcat(values(Aadj_m2p_remaining)...))) - length(Pdiscarded)

@test poles_used_as_per_Am2p == poles_used_as_per_Am2p_remaining_and_Pdiscarded == poles_used_as_per_Aadj_m2p_remaining_and_Pdiscarded == poles_used

@unpack preprocess1_steps, preprocess2_steps, preprocess3_steps = graph;
HF.myprintln(true, "preprocess1_steps=$(preprocess1_steps)")
HF.myprintln(true, "preprocess2_steps=$(preprocess2_steps)")


# open("profile_summary.txt", "w") do io
#     Profile.print(io; format=:flat, sortedby=:count)
# end