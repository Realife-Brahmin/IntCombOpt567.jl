# main.jl
include("./src/setupIntCombOpt567.jl")
Revise.track(IntCombOpt567.setCoveringHeuristics)

matDict = sCH.txt2mats("rawData/project02/p4m6.txt")