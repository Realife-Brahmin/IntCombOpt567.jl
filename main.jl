# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

# testCase = "p4m6"
# testCase = "phase1"
testCase = "cap360"

# sim_time = @belapsed begin
    g = sCH.initializeGraph("rawData/project02/" * testCase * ".txt")
    @profile sCH.solveSetCoveringProblem!(g)
    global graph = g
# end

@unpack poles_used, meters_covered, Acov, m = graph
@test meters_covered == m
myprintln(true, "value=$(poles_used), time=$(sim_time)")

# display(Acov)