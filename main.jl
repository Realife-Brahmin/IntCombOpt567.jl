# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

testCase = "p4m6"
# testCase = "phase1"
# testCase = "cap360"

graph = sCH.initializeGraph("rawData/project02/" * testCase * ".txt")

graph = sCH.solveSetCoveringProblem(graph)

@unpack poles_used, meters_covered, Acov, m = graph;
@test meters_covered == m
myprintln(true, "Before value=$(poles_used)")
display(Acov)

sCH.removePole!(graph, 2, verbose=true)
@unpack meters_covered, poles_used, Acov, m = graph;
@test meters_covered == m 
myprintln(true, "After value=$(poles_used)")
display(Acov)