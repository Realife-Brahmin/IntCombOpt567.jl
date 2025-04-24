module IntCombOpt567

export setCoveringHeuristics

using AmplNLWriter
using Gurobi
using HiGHS
using JuMP
using Parameters: @unpack, @pack!
using SparseArrays

include("./setCoveringHeuristics.jl")

using .setCoveringHeuristics

export txt2mats

end
