module IntCombOpt567

export helperFunctions
export setCoveringHeuristics

using AmplNLWriter
using Gurobi
using HiGHS
using JuMP
using Parameters: @unpack, @pack!
using SparseArrays

include("./helperFunctions.jl")
include("./setCoveringHeuristics.jl")

using .helperFunctions
using .setCoveringHeuristics

export
    chooseNextPole,
    getNextPole,
    myprintln,
    solveSetCoveringProblem,
    txt2mats

end
