push!(LOAD_PATH, "./")
push!(LOAD_PATH, "./src/")

using SLgraphene
using Constants
using Printf

function mkfolder(path)
	if(isdir(path))
		rm(path, recursive=true)
	else
		mkdir(path)
	end
end




λ = 84; σ = 6.5*nm; Δh0 = 2.5*nm; fᵤ = 1; U = 0*eV; μ = 0*eV; V₀ = 0.0*eV; U₂ = 0.00*eV
#main(λ, σ, Δh0, fᵤ, U, μ, U₂)

toppath = "./testing/10-17-21-bigruns"
save = true
rm(toppath, recursive=true)
mkdir(toppath)
#σ1 = 0.5*nm; σ2 = 3.5*nm; dσ = 1*nm
#h0 = 1*nm; h2 = 5*nm; dh = 1*nm
for σ = 4*nm:1*nm:7*nm
	path = toppath*"/"*(@sprintf("L#%i_S#%g_dh0#%g_fu#%g_U#%g_mu#%g_V0#%g_U2#%g",λ,σ,Δh0,fᵤ,U,μ,V₀,U₂))
	mkdir(path)
	#for σ = σ1:1*nm:σ2
	main(λ, σ, Δh0, fᵤ, U, V₀, V₀, U₂,save,path)
	#end
end


