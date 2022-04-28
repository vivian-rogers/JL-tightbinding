push!(LOAD_PATH, "./src/")
push!(LOAD_PATH, "../src/")

module Driver
#using Plots
#using PlotStuff
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
#using Nanospheres
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
function runBands(p,nk, H, Q, proj::Bool=false,save=false,path="",name="")
	
	#default list of interesting points -- notation can be changed if desired, look in 
	klist = ["Γ","M","X","Γ","Z","A","R","Z"]
	#klist = ["Γ","M","X","Γ","A","M","X","R","X","Γ","Z"]
	#klist = ["Γ","M","X","Γ","Z","A","R","Z"]

	# generates k-name -> k-value correspondence 
	kdict = p.kdict

	println("\n========= Entering Bands calculation =========")
	println("Getting eigenvalues of between k = ")
	show(klist)
	println("...")
	
	# diagonalizes hamiltonian on k-path in middle of spectrum, see getbands for # of eigvls
	E, Estates = getBands(klist, kdict, nk, p.a, H)
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

# takes in the parameters defined in runs.jl, returns H(k)
function runSCF(p,A,save=false,path="",name="")
	
	# returns static part of hamiltonian, periodic boundary hopping structs, up/down spin density, 
	# position values of sites, pseudo B field(R)
	# (mess with spin polarization in /src/ConstructHamiltonian.jl, needs testing. Currently n₊ = n₋)
	H = Hgen(p,A)
	#H₀, edgeNN_arr, n₊, n₋, R, pseudoB = Hgen(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,θ,save,path)
	
	# older implementation -- may be more efficient to gen H(k) in Hgen
	# H, n₊, n₋, R, pseudoB = Hgen(λ,Δh0,σ,fᵤ,U,μ,V₀,U₂,θ,save,path)
	#=
	# generates stats on the pseudo B field
	Brms, Bave, Bmax = pseudoBstats(pseudoB)
	println("pseudo-B stats: ")
	@printf("<|B|> = %g T, √(B²) = %g T, max(B) = %g T\n", Bave, Brms, Bmax)
	

	
	# prepping for final H definition
	# this is nontrivial to read, look into map function in functional languages for explanation
	edgeNNs = [NN for NN in edgeNN_arr]
	
	# the grid-site values of the relevant periodic edge bonds
	aVals = map(NN -> Int(NN.a),edgeNNs)
	bVals = map(NN -> Int(NN.b),edgeNNs)
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	function H(k)
		
		# connect the edges using bloch's theorem
		# needs phase offset if breaking time reversal symmetry, or perhaps define in Hgen
		tVals = map(NN->ComplexF64(NN.t*exp(im*k⋅(NN.N[1]*a₁+NN.N[2]*a₂))), edgeNNs)
		
		# generates the matrix for cₐ†cᵦ for bonds leading out of SL primitive unit cell
		H_edge = sparse(aVals,bVals,tVals)

		# returns total hamiltonian
		return H₀ .+ H_edge
	end

	# can plot the spin densities !
	if(save)
		fig₊ = plotScatter(R,n₊,"x position (nm)","y position (nm)","|ψ|²");
		fig₋ = plotScatter(R,n₋,"x position (nm)","y position (nm)","|ψ|²");
		SaveFigure(fig₊,path,"ψ²_u"*name)
		SaveFigure(fig₋,path,"ψ²_d"*name) 
	end=#
	println("returning from SCF")
	# returns call to H(k), Brms value for plotting and stats
	return H
end

function A(R::Vector{Float64})
	return [0;0;0]
end

# What Is To Be Done?
function main(p,save=false,path="")
	
	# number of sites
	println("=========== InBi ==============")
	println("Running simulation with parameters: ")
	println("params = $p") # oh this is so ugly, clean it later
	if(save) println("Saving outputs to $path\n") end

	println("Generating static hamiltonian...")
	H = runSCF(p,A,save,path)

	# some arbitrary operators for bands routine
	
	# project onto Q = |In><In| for bands purposes
	# [unit cell] ⊗ [A/B site] ⊗ [atom type] ⊗ [px, py] ⊗ [spin]
	Q = I(p.n)⊗I(2)⊗(Diagonal([1 0]))⊗I(p.norb)⊗I(2)
	runBands(p,2^8,H,Q,true)
	#DOS, Evals = runDOS(20,H,λ,save,path,Beff)
	#runLDOS(20, H, λ,save,path,true,Beff)
	println("Done!\n")

	# if sweeping, can return something here
end
end
