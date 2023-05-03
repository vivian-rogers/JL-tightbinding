

module DOS
using UsefulFunctions
#using Arpack
using LinearAlgebra
using Arpack
using GMaterialParameters
using Logging
using PlotlyJS
using Distributed
using ProgressBars
using StatsBase
using Operators
#using GLMakie

export energySurface, kslice, getDOS, getLDOS, eigSurface


function genBZmesh(B::Matrix,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    X1 = B*[1/2;0;0];
    X2 = B*[0;1/2;0];
    X3 = B*[0;0;1/2];
    kx =  LinRange(-X1[1],X1[1],2*nx+1)
    ky =  LinRange(-X2[2],X2[2],2*ny+1)
    kz =  LinRange(-X3[3],X3[3],2*nz+1)
    return kx, ky, kz
end

function wilsonLoop(p::NamedTuple, H::Function)


function nonAbelianBerryCurvature(p::NamedTuple, H::Function, k::Vector{Float64}, dk::Vector{Float64})
        
end
