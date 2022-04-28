push!(LOAD_PATH, "./src/")

module system
using Plots
using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using MaterialParameters
using Printf
using FastClosures
using DOS
using ConstructHamiltonian
using Bands

export main

function runBands(nk,H,Q,proj=false,save=false)
	klist = ["Γ", "X", "W", "K", "Γ", "L", "W","U","X"]
	println("\n========= Entering Bands calculation =========")
	println("Getting eigenvalues of between k = ")
	show(klist)
	println("...")
	E, Estates = getBands(klist, kdict, nk, a, H)
	println("Plotting...")
	if(proj)
		projStates = expectedValue(Q,Estates)
		fig = plotBands(klist,nk,E, projStates)
	else
		fig = plotBands(klist,nk,E,2)
	end
	if(save) SaveFigure(fig,path,"bands"*name) end
end

function runLDOS(n, H, λ,save=false,path="",name="")
	klist, weights = kgrid(n,λ) 
	println("\n========= Entering LDOS calculation =========")
	println("Binning eigvls over γ-κ-κ' IBZ, density = $n")
	LDOS, DOS = getLDOS(klist, weights, H, λ)
	println("Plotting...")
	fig1 = plot2D(LDOS,"E (eV)", "DOS")
	#fig2 = plot1D(Evals,DOS,"E (eV)", "DOS")
	if(save) SaveFigure(fig1,path,"LDOS"*name) end
	#if(save) SaveFigure(fig2,path,"DOS"*name) end
end

function runDOS(n, H, λ,save=false,path="", Bavg=0)
	klist, weights = kgrid2(n,λ) 
	println("\n========= Entering DOS calculation =========")
	println("Binning eigvls over γ-κ-κ' IBZ, density = $n")
	DOS, Evals = getDOS(klist, weights, H)
	println("Plotting...")
	fig = plot1D(Evals,DOS,"E (eV)", "DOS")
	if(Bavg > 0)
		fig = addLLs(Bavg, fig, 5, 0, false, true)
	end
	if(save) SaveFigure(fig,path,"DOS") end
	return DOS, Evals
end

#multiplicity of a_SL/a_graphene, stddev of NS, Δheight of NS profile
#λ = 30; σ = 2*nm; Δh0 = 2*nm; fᵤ = 1; U = 0.3*eV; μ = 0.01*eV
#N = 2*λ^2 #total number of atoms in system
#projKet = I(Int(N/2))⊗[0;1] #project onto B site

#x, y, z = NSprofileGen(λ,Δh0,σ)
#surface(x,y,z)

function runSCF(spin,save=false,path="")
	#H, n₊, n₋, R, pseudoB = Hgen(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,save,path)
	#H₀, edgeNNs, n₊, n₋, R, pseudoB = Hgen(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,θ,save,path)
	#H, n₊, n₋, R, pseudoB = Hgen(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,θ,save,path)

	H = Hgen(spin)
	#println("Generating H(k) from H_static and nearest neighbors")
	#=function H(k::Array{Int64,1}) 		
		let
			Hstatic = Hstatic
			edgeNNs = edgeNNs
			λ = λ
		end
		return makeH(Hstatic,edgeNNs,λ,k)
	end =#
	#H(k::Array{Int64,1}) = @closure k -> makeH(H₀,edgeNNs, λ,k)
	#function H(k)
	#	return makeH(H₀, edgeNNs, λ, k)
	#end
	#H(k) = @closure k -> (H₀ .+ sparse(map(NN->NN.a,edgeNNs), in edgeNNs)
	#H(k) = @closure k -> (H₀ .+ sparse(map(NN->NN.a,edgeNNs),map(NN->NN.b,edgeNNs),map(NN->NN., in edgeNNs)
	#show(H([0;0;0]))
	#show(makeH(k,H₀, edgeNNs, λ))
	#show(H([0;0]))
	return H
	#return Hstatic, edgeNNs, Brms
end


function main(spin=true,save=false,path="")
	N = 4
	println("=========== Heusler ==============")
	println("Running simulation with parameters: ")
	println("a = $(a*10^10)*Å")
	if(save) println("Saving outputs to $path\n") end

	println("Generating static hamiltonian...")
	#H, Brms = runSCF(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,θ,save,path)
	H = runSCF(spin,save,path)
	#H(k) = k -> makeH(Hstatic,edgeNNs,λ,k)
	#Q = closestPeak(λ) 
	if(spin==true)
		Q = I(N)⊗σ₃ # spin operator 
	else
		Q = I(2)⊗[1 0 ; 0 0]
	end
	#Q = I(N)⊗σ₃ # spin operator 
	#Q = I(N)⊗σ₃ # spin operator 
	#runBands(2^7,H,false)
	#show(size(H([0;0])))
	#show(size(Q))
	runBands(2^8,H,Q,true,false)
	#DOS, Evals = runDOS(20,H,λ,save,path,Brms)

	println("Done!\n")
end

main(true)

end
