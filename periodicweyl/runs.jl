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

#nx = 10; ny = 10; nz = 10; 
nx = 200; ny = 1; nz = 1; 
#nx = 12; ny = 1; nz = 1; 
# superlattice basis vectors, in basis of a_1, a_2, a_3
SL1 = [nx; 0; 0]; SL2 = [0; ny; 0]; SL3 = [0; 0; nz]

runparams = (
             μ = 0.0*eV, μ_disorder = 0.025*eV, 
             #E_samples = [0.1,0.2,0.3,0.4],
             E_samples = [0.05],
             #E_samples = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
             #E_samples = [E for E = 0:0.1:1.5],
             #E_samples = [E for E = -6:0.5:6],
             #E_samples = [E for E = -6:0.005:6],
             nk = 80,
             β = 0.10*eV, runtype = "neellattice", fieldtype = "β", η = 10^-4, ηD = 10^-3, 
             verbose = false, plotfield = true, bands=false, θ=350.0, sweep="none",
             electrodeMagnetization=true,electrodeMaterial="weyl",
             #electrodeMagnetization=false,electrodeMaterial="metal",
             deviceMagnetization=true,deviceMaterial="weyl",
             startDWs = 20*nm, DWwidth = 18*nm, DWspacing = 100*nm, 
             #startDWs = 30*nm, DWwidth = 3*nm, DWspacing = 10*nm, 
             λ = 2*nm,
             #electrodeMagnetization=true,electrodeMaterial="mtjweyl",
             #deviceMagnetization=false,deviceMaterial="ins",
             prune=[]
             #prune=["z"]
            )

parms = merge(params,runparams)
p = genSL(parms, nx, ny, nz, SL1, SL3, SL3, parms.runtype, parms.fieldtype) # generate SL params


function runFieldTexture(p::NamedTuple)
    if(p.fieldtype=="A")
            A = Agen(p,p.runtype,10^8*4*μₑ)
    else
            A = βgen(p,p.runtype,p.β,p.θ,p.startDWs)
    end
    return main(p,A)
end

function nxtoArg(nx::Int)
    SL1 = [nx; 0; 0]
    return genSL(parms, nx, ny, nz, SL1, SL3, SL3, parms.runtype, parms.fieldtype) # generate SL params
end

function startDWstoArg(startDWs::Float64)
    return merge(p, (sweep = "T(Stripe Domain Position,E)", startDWs = startDWs))
end

function θtoArg(θ::Float64)
    return merge(p, (sweep = "T(DW angle,E)", θ = θ))
end

runFieldTexture(p)

#Sweep1DSurf(runFieldTexture,θtoArg,[θ for θ = 0:10.0:180],p.E_samples,"θ DW angle (degrees)", "Energy (eV)", "T (e²/h)")
#Sweep1DSurf(runFieldTexture,startDWstoArg,[DWstart for DWstart = -5*nm:(1*nm):(0.7*p.SLa₁[1])],p.E_samples,"Stripe DW start position (nm)", "Energy (eV)","T (e²/h)",(1/nm),false)
#Sweep1DSurf((T -> log10.(T))∘runFieldTexture,startDWstoArg,[DWstart for DWstart = -5*nm:(5*nm):(30*nm)],p.E_samples,"Stripe DW start position (nm)", "Energy (eV)","T (e²/h)",(1/nm),false)
#Sweep1DSurf((T -> log10.(T))∘runFieldTexture,θtoArg,[θ for θ = 0.0:10.0:180.0],p.E_samples,"θ DW angle (degrees)", "Energy (eV)","log₁₀(T) (e²/h)",false)

end
