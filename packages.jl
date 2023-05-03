using Pkg

dependencies = [
    "IJulia",
    "LinearAlgebra",
    "Printf",
    "SparseArrays",
    "SparseSuite",
    "Plots",
    "Makie",
    "JLD2",
    "Symbolics",
    "BlockDiagonals",
    "Distributions",
    "Arpack",
    "Dates",
    "Optim",
    "Distributed",
    "LaTeXStrings",
    "PyCall",
    "WebIO",
    "ForwardDiff",
    "StatsBase",
    "ProgressBars",
    "Logging",
    "DelimitedFiles",
    "Random",
    "Makie",
    "GLMakie",
    "ColorSchemes",
    "Suppressor",
    "Measures"
]

Pkg.add(dependencies)
