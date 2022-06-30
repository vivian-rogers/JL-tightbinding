push!(LOAD_PATH, "./src/")
push!(LOAD_PATH, "../src/")

module Driver
#using Plots
#using PlotStuff
using InBi
using Electrodes
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
#using GMaterialParameters
using Printf
#using FastClosures
using DOS
using SparseArrays
using ConstructHamiltonian
using InBiNNs
using VectorPotential
#using Nanospheres
using NEGF
using Bands
using PlotStuff

export main


# runs the bands subroutine, see /src/bands.jl. Sweeps along k-path, gets eigenvalues, (optional) projects eigenstates, returns figure
# input: number of interpolated k-points between high-sym points, H(k), λ, operator Q (same size as H),
# bool (project onto Q? y or n), bool (save bands to file?), path, name
# 
# output: figure, with bands projected (or not))
# can project eigstates onto an operator Q that is defined in the hilbert space of H(k) at a defined k.
# Needs a new subroutine to calculate berry curvature or chern number  
function runBands(p,nk, H, Q, proj::Bool=false,arpack::Bool=false,save=false,path="",name="")
	
	#default list of interesting points -- notation can be changed if desired, look in 
	#klist = ["Γ","M","X","Γ","Z","A","R","Z","A","M","X","R","X","Γ","Z"]
	#klist = ["Γ","M","X","Γ","A","M","X","R","X","Γ","Z"]
	#klist = ["X₁","Γ","Z"]
	#klist = ["X","Γ","-X"]
	#klist = ["X₁","Γ","-X₁"]
	#klist = ["Γ","M","X₁","Γ","-X₁"]
	#klist = ["Γ","M","X₁","Γ","-X₁","X₂","Γ","-X₂"]
	#klist = ["Γ","M","X","Γ","-X","Γ","Z","A","R","Z"]

	# generates k-name -> k-value correspondence 
	klist = p.klist
	kdict = p.kdict

	println("\n========= Entering Bands calculation =========")
	println("Getting eigenvalues of between k = ")
	show(klist)
	println("...")
	
	# diagonalizes hamiltonian on k-path in middle of spectrum, see getbands for # of eigvls
	E, Estates = getBands(klist, kdict, nk, p.a, H, arpack)
	println("Plotting...")
	#Q = closestPeak(λ) 
	#Q = I(2)⊗(Diagonal([1 0]))⊗I(2)⊗I(2)
	#Q = I(2)⊗(Diagonal([1 0]))⊗Diagonal([1 1])⊗I(2)
	if(proj)
		# projects onto Q
		projStates = expectedValue(Q,Estates)
		fig = plotBands(klist,nk,E, projStates)
	else
		fig = plotBands(klist,nk,E,2)
	end
	if(save) SaveFigure(fig,path,"bands"*name) end
end


# gets LDOS by projecting eigenstates of hamiltonian onto sites, over k-grid. See /src/DOS.jl
# input: n = density of evenly-spaced k-grid (0 = γ, 1 = γ, κ', κ, 2 = etc) -- 20 gives reasonable results
# H(k), λ, boolean to save figs or not, path, boolean to plot total DOS, pseudo-B field to plot pLLs
# Not super efficient, needs rewriting in recursive green's function formalism in getLDOS
function runLDOS(n, H, λ,save=false,path="",plotDOS=true,Bavg=0)
	# generates the evenly-spaced k-grid to integrate spectrum over
	# kgrid2 includes redundant k-points but IBZ grid generator was bugged, needs fixing for optimization.
	klist, weights = kgrid2(n,λ) 
	println("\n========= Entering LDOS calculation =========")
	println("Binning eigvls over γ-κ-κ' IBZ, density = $n")
	LDOS, DOS, Evals = getLDOS(klist, weights, H, λ)
	println("Plotting...")
	
	# generate figure for path along diagonal of superlattice primitive unit cell
	fig1 = plot2D(LDOS,0,nm^-1*√(3)*λ*a,minimum(Evals),maximum(Evals),"Distance along (a₁+a₂) (nm)","E (eV)")
	
	if(save) SaveFigure(fig1,path,"LDOS") end
	if(plotDOS) 
		println("Plotting...")
		fig = plot1D(Evals,DOS,"E (eV)", "DOS")
		if(Bavg > 0)
			fig = addLLs(Bavg, fig, 5, 0, false, true)
		end
		if(save) SaveFigure(fig,path,"DOS") end
	end
end

# similar to LDOS above, less functionality 
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

# function to calculate transport from device parameters
function NEGF_2contacts_1layer(p::NamedTuple,A::Function)
        println("============= NEGF transport ============")
        Electrodes = [
		Electrode([-1,0],[0,p.ny],[0,p.nz],p.ny*p.nz,"-x","weyl",A)
		Electrode([p.nx,p.nx+1],[0,p.ny],[0,p.nz],p.ny*p.nz,"+x","weyl",A)
	]
	negf_params = (prune = union(["x"],(p.prune)), nelectrodes=2)
	#negf_params = (prune = union(["x"],(p.prune)), nelectrodes=2, η=10^-8*eV)
	negf_params = merge(p,negf_params)
	H = runSCF(negf_params,A) # generates H(k)
	println("Generating self-energy matrices for electrodes...")
	Σks = genΣₖs(p,Electrodes) #i.e. the Σᵢ(E) at a given k value. Call with Σₖ = Σks(k); then Σ = Σₖ(E)
	println("Defining Gʳ, Current operator, etc...")
	genGʳ, genT, genA = NEGF_prep(negf_params,H,Σks) # returns the functions to generate [quantity](E) by calling genQ(k)
	# let's sample at
	nk = 60
	kgrid, kweights, kindices = genBZ(negf_params,0,nk,nk) # generate surface BZ points
	E_samples = [0.5]
	#E_samples = [0.2,0.5,1.0,2.0,5.0]
	#E_samples = [0.5]
	#E_samples = [E for E = 0.1:0.5:2.1]
	#E_samples = [E for E = 0.1:1.0:2.1]
        TofE, Tmap, TmapList = totalT(genT, kindices, kgrid, kweights, E_samples, minimum(E_samples))
        #display(TmapList)
	#display(Tmap)
	display(TofE)
	plot1D(TofE,E_samples,"Transmission","E (eV)",0.0,10.0,minimum(E_samples),maximum(E_samples))
        plotMat(Tmap',"k₂","k₃")
        #plotMat([i[1] for i in kindices],[i[2] for i in kindices],TmapList,"k₂","k₃")
        #plot2D([i[1] for i in kindices],[i[2] for i in kindices],TmapList,"k₂","k₃")
	#plotSurf([k[1] for k in kgrid], [k[2] for k in kgrid], Tmap, "ky (π/a₂)", "kz (π/a₃)")	
	#plotSurf([k[1] for k in kgrid], [k[2] for k in kgrid], Tmap, "ky (π/a₂)", "kz (π/a₃)")	
	

end

# takes in the parameters defined in runs.jl, returns H(k)
function runSCF(p,A,save=false,path="",name="")
	H = Hgen(p,A)
	println("returning from SCF")
	# returns call to H(k), Brms value for plotting and stats
	return H
end

function A(R::Vector{Float64})
	return [0;0;0]
end

# What Is To Be Done?
function main(p,A=A,save=false,path="")
	
	# number of sites
	println("=========== InBi ==============")
	println("Running simulation with parameters: ")
	println("params = $p") # oh this is so ugly, clean it later
	if(save) println("Saving outputs to $path\n") end

	println("\n")
	println("Generating static hamiltonian...")
	H = runSCF(p,A,save,path)

	testH = H([0;0;0])
	println("size of H = $(size(testH))")
	#println("testH₀ = \n $testH")
	# some arbitrary operators for bands routine
	
	# project onto Q = |In><In| for bands purposes
	# [unit cell] ⊗ [A/B site] ⊗ [atom type] ⊗ [px, py] ⊗ [spin]
	Q = I(p.n)⊗I(p.nsite)⊗I(p.norb)⊗σ₂
        #Q = distFromDW(p,RvalsGen(p))⊗I(p.norb)⊗(2) 
        #Q = zpos(RvalsGen(p))⊗I(p.norb)⊗I(2)
	#runBands(p,2^6,H,Q,true,p.arpack)
	NEGF_2contacts_1layer(p,A)
	#DOS, Evals = runDOS(20,H,λ,save,path,Beff)
	#runLDOS(20, H, λ,save,path,true,Beff)
	println("Done!\n")

	# if sweeping, can return something here
end
end
