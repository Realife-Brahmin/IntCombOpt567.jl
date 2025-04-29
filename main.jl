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
# scoring_function = "greedy"
# scoring_function = "score1"
scoring_function = "score2"
benchmarkTime = false
# benchmarkTime = true

g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats,
    scoring_function=scoring_function)

#region solveSetCoveringProblem
if benchmarkTime
    # just time the solve call
    sim_time = @belapsed begin
        g_local = deepcopy(g)
        sCH.solveSetCoveringProblem!(g_local)
    end
    println("Solved in $(sim_time) s")
else
    sim_time = "not_benchmarked"
end

# profile separately (only once)
@profile begin
    sCH.solveSetCoveringProblem!(g)
    global graph = g
end
#endregion solveSetCoveringProblem

@unpack poles_used, meters_covered, Acov, m = graph
@test meters_covered == m
myprintln(true, "value=$(poles_used)")

# display(Acov)

open("profile_summary.txt", "w") do io
    Profile.print(io; format=:flat, sortedby=:count)
end