# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

testCase = "p4m6"
# testCase = "phase1"
# testCase = "cap360"

graph = sCH.initializeGraph("rawData/project02/" * testCase * ".txt")

graph = sCH.solveSetCoveringProblem(graph)

@unpack poles_used, Acov = graph;
myprintln(true, "Before value=$(poles_used)")
display(Acov)
graph = sCH.removePole(graph, 2, verbose=true)


@unpack poles_used, Acov = graph;
myprintln(true, "After value=$(poles_used)")
display(Acov)