# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)
Revise.track(IntCombOpt567.helperFunctions)

testCase = "p4m6"
# testCase = "phase1"
# testCase = "cap360"

graphState = sCH.initializeGraph("rawData/project02/" * testCase * ".txt")

graphState = sCH.solveSetCoveringProblem(graphState)

@unpack poles_used, Acov = graphState;
myprintln(true, "Before value=$(poles_used)")
display(Acov)
graphState = sCH.removePole(graphState, 2, verbose=true)


@unpack poles_used, Acov = graphState;
myprintln(true, "After value=$(poles_used)")
display(Acov)