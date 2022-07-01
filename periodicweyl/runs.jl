push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./periodicweyl/")
push!(LOAD_PATH, "./src/")


module runs
using InBi
using SLutils
using Driver
using Constants
using UsefulFunctions
using LinearAlgebra
using PlotStuff
#main(params)

nx = 10; ny = 1; nz = 1; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]

runparams = (
             μ = 0.0*eV, μ_disorder = 0.0*eV, 
             #E_samples = [0.1,0.2,0.3,0.4],
             E_samples = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
             nk = 80,
             β = 0.5*eV, runtype = "blochdw", fieldtype = "β", η = 3*10^-3, 
             verbose = false, plotfield = false, bands=false, θ=30.0, sweep="none"
            )

parms = merge(params,runparams)
p = genSL(parms, nx, ny, nz, SL1, SL3, SL3, parms.runtype, parms.fieldtype) # generate SL params


function runFieldTexture(p::NamedTuple)
    if(p.fieldtype=="A")
            A = Agen(p,p.runtype,10^8*4*μₑ)
    else
            A = βgen(p,p.runtype,p.β,p.θ)
    end
    return main(p,A)
end

function θtoArg(θ::Float64)
    return merge(p, (sweep = "T(DW angle,E)", θ = θ))
end

#runFieldTexture(p)

Sweep1DSurf(runFieldTexture,θtoArg,[θ for θ = 0.0:10.0:90.0],p.E_samples,"θ DW angle (degrees)", "Energy (eV)")

end
