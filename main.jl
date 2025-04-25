# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

# testCase = "p4m6"
testCase = "phase1"
# testCase = "cap360"
cleanupRepeats = 1
# cleanupRepeats = 10

sim_time = @belapsed begin
    g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt",
    cleanupRepeats=cleanupRepeats)
    @profile sCH.solveSetCoveringProblem!(g)
    global graph = g
end 
myprintln(true, "sim_time=$(sim_time)")

@unpack poles_used, meters_covered, Acov, m = graph
@test meters_covered == m
myprintln(true, "value=$(poles_used)")

# display(Acov)

# Save the top 20 lines (most sampled = most time spent) to a file
open("profile_summary.txt", "w") do io
    Profile.print(io; format=:flat, sortedby=:count)
end