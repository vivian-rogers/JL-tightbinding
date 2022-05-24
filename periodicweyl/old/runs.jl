push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

using SLgraphene
using Constants
using Printf
using PlotStuff

function mkfolder(path)
	exists = isdir(path)
	if(exists)
		println("$path already exists...")
		#rm(path, recursive=true)
	else
		mkdir(path)
	end
	return exists
end




const λ = 10; const σ = 0.5*nm; const Δh0 = 0.1*nm; 
const θ = 0; fᵤ = 1; const U = 0*eV; const μ = 0*eV; 
const V₀ = 0.00*eV; const U₂ = 0.00*eV
#main(λ, σ, Δh0, fᵤ, U, μ, U₂)

const toppath = "./testing/11-5-2021"
const save = false
#rm(toppath, recursive=true)
mkfolder(toppath)
#σ1 = 0.5*nm; σ2 = 3.5*nm; dσ = 1*nm
#h0 = 1*m; h2 = 5*nm; dh = 1*nm
const hVals = [h*nm for h = 0:0.5:1]
sweepDOS = zeros(size(hVals)[1],100)
for i in eachindex(hVals)
	Δh0 = hVals[i]
	path = toppath*"/"*(@sprintf("L-%i_S-%g_dh0-%g_fu-%g_theta-%g_U-%g_mu-%g_V0-%g_U2-%g",λ,σ,Δh0,fᵤ,θ,U,μ,V₀,U₂))
	mkfolder(path)
	#if(mkfolder(path)==true)
	#for σ = σ1:1*nm:σ2
	#main(λ, σ, Δh0, fᵤ, U, μ, V₀, U₂, θ, save,path)
	DOS, Evals = main(λ, σ, Δh0, fᵤ, U, μ, V₀, U₂, θ, save,path)
	sweepDOS[i,:] .= DOS
	#end
end

fig = plot2D(sweepDOS')


