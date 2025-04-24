# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)

# testCase = "p4m6"
testCase = "phase1"
# testCase = "cap360"

matDict = sCH.txt2mats("rawData/project02/" * testCase * ".txt")
