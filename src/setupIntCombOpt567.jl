using BenchmarkTools
using Revise
using IntCombOpt567
using Parameters
using Profile
using ProfileView
using SparseArrays
using Test

import .IntCombOpt567.setCoveringHeuristics as sCH
import .IntCombOpt567.helperFunctions as HF

Revise.revise()