# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

testCase = "p4m6"
# testCase = "phase1"
# testCase = "cap360"

graphState = sCH.txt2graph("rawData/project02/" * testCase * ".txt")

graphState = sCH.solveSetCoveringProblem(graphState)

@unpack poles_used = graphState;
myprintln(true, "value=$(poles_used)")