# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

# testCase = "p4m6"
# testCase = "phase1"
testCase = "cap360"
cleanupRepeats = 1
# cleanupRepeats = 10
scoring_function = "greedy"
# scoring_function = "score1"
# scoring_function = "score2"
benchmarkTime = false
# benchmarkTime = true

g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats,
    scoring_function=scoring_function)

#region solveSetCoveringProblem
if benchmarkTime
    sim_time = @belapsed begin
        @profile sCH.solveSetCoveringProblem!(g)
        global graph = g
    end
else
    @profile begin
        global graph = sCH.solveSetCoveringProblem!(g)
    end
    sim_time = "not_benchmarked"
end 
#endregion solveSetCoveringProblem

myprintln(true, "sim_time=$(sim_time)")

@unpack poles_used, meters_covered, Acov, m = graph
@test meters_covered == m
myprintln(true, "value=$(poles_used)")

# display(Acov)

# Save the top 20 lines (most sampled = most time spent) to a file
open("profile_summary.txt", "w") do io
    Profile.print(io; format=:flat, sortedby=:count)
end