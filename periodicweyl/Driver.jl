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
using MakiePlots
using SaveStuff
using JLD2

export main


# runs the bands subroutine, see /src/bands.jl. Sweeps along k-path, gets eigenvalues, (optional) projects eigenstates, returns figure
# input: number of interpolated k-points between high-sym points, H(k), λ, operator Q (same size as H),
# bool (project onto Q? y or n), bool (save bands to file?), path, name
# 
# output: figure, with bands projected (or not))
# can project eigstates onto an operator Q that is defined in the hilbert space of H(k) at a defined k.
# Needs a new subroutine to calculate berry curvature or chern number  
function runBands(p,nk, H, Q, proj::Bool=false,arpack=false,save=false,path="",name="")
	
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
	E, Estates = getBands(klist, kdict, nk, p.a, H, 32)
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
        if(p.savedata) 
            #jldsave(p.path * "bands.jld2",fig) 
            mkdelim(p.path * "bandsdata.txt", [collect(1:size(E)[1]) E projStates])
        end
        if(p.save) SaveFigure(fig,p.path,"bands"*name) end
end


# gets LDOS by projecting eigenstates of hamiltonian onto sites, over k-grid. See /src/DOS.jl
# input: n = density of evenly-spaced k-grid (0 = γ, 1 = γ, κ', κ, 2 = etc) -- 20 gives reasonable results
# H(k), λ, boolean to save figs or not, path, boolean to plot total DOS, pseudo-B field to plot pLLs
# Not super efficient, needs rewriting in recursive green's function formalism in getLDOS
#=

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
	return DOS, Eval
end
=#

# function to calculate transport from device parameters
function NEGF_2contacts_1layer(p::NamedTuple,A::Function,returnvals)
        println("============= NEGF transport ============")
        Electrodes = [
		Electrode([-1,0],[0,p.ny],[0,p.nz],p.ny*p.nz,"-x",p.electrodeMaterial,A);
		Electrode([p.nx,p.nx+1],[0,p.ny],[0,p.nz],p.ny*p.nz,"+x",p.electrodeMaterial,A)
	]
        plot_params = (prune = setdiff(p.prune,["x"]), sweep="plot"); 
        plot_params = merge(p,plot_params);
        display(plot_params)
        negf_params = (prune = union(["x"],(p.prune)),verbose=false, nelectrodes=size(Electrodes)[1])
	negf_params = merge(p,negf_params)
	H = runSCF(negf_params,A,returnvals) # generates H(k)
	println("Generating self-energy matrices for electrodes...")
	Σks = genΣₖs(plot_params,Electrodes) #i.e. the Σᵢ(E) at a given k value. Call with Σₖ = Σks(k); then Σ = Σₖ(E)
	println("Defining Gʳ, Current operator, etc...")
	genGʳ, genT, genA, genScatteredT = NEGF_prep(negf_params,H,Σks) # returns the functions to generate [quantity](E) by calling genQ(k)
	# let's sample at
        γ⁵ = I(p.nx*p.ny*p.nz)⊗τ₁⊗σ₀;
        if(p.mixedDOS==true)
            mdE = p.E_samples[1]*eV; η = 10^-(3.0)
            function plottingGʳ(k::Vector{Float64})
                function Gʳ(E::Float64)
                    Σₗ = Σks[1]; Σᵣ = Σks[2]
                    return pGrInv((E+im*η)*I(p.nx*p.ny*p.nz*p.norb*2) .- H(k) .- Σₗ(k)(E) .- Σᵣ(k)(E),4,false) 
                end
                return Gʳ
            end
            # left and right-handed states
            Operators = [(1/2)*(I(p.nx*p.ny*p.nz*p.norb*2).+d*γ⁵)  for d = [-1,+1]]
            DOS = sitePDOS(plot_params,plottingGʳ,Operators, mdE)
			#G = genGʳ([0.1;0.1;0.1])(0.02*eV)
			#Test = DOS([0.1;0.1;0.1]);
            #testDOS(k) = ones(p.nx)*(k⋅k)/nm^2;
            nkDOS = 250
            fsfig = mixedDOS(plot_params,DOS,nkDOS,nkDOS)
            if("mixedDOS" ∈ p.returnvals)
                push!(returnvals,fsfig)
            end
        end
        nkx = p.nk*!("x" ∈ negf_params.prune); nky = p.nk*!("y" ∈ negf_params.prune); nkz = p.nk*!("z" ∈ negf_params.prune);
	
        kgrid, kweights, kindices, kxs, kys, kzs = genBZ(negf_params,nkx,nky,nkz) # generate surface BZ points
        println("Sweeping transmission over kgrid: $(nkx*2+1), $(nky*2+1), $(nkz*2+1) ...")
        #TofE, Tmap, TmapList = totalT(genT, kindices, 0.3 .* kgrid, kweights, p.E_samples, minimum(p.E_samples))
        #TofE, Tmap = totalT(genT, kindices, 0.3 .* kgrid, kweights, p.E_samples, minimum(p.E_samples))
        parallelk = ((nkx+1)*(nky+1)*(nkz+1) > 8)
        if(p.nk > 0)
            S = 0.15 # scale for k-map
        else
            S = 1
        end
        #println("parallelk = $parallelk, negf_params.prune = $(negf_params.prune)")
        Operators = [I(p.nx*p.ny*p.nz*p.norb*2), γ⁵]
        TofE, Tmap = totalT(genScatteredT, kindices, S .* kgrid, kweights, p.E_samples, p.E_samples[1], parallelk, Operators)
        TofE = S^2*TofE
        figh = pyplotHeatmap(S*kys/(π/p.a),S*kzs/(π/p.a),Tmap',"ky (π/a)","kz (π/a)","T(ky,kz)",:nipy_spectral, p.savedata, p.path)
        if("tplot" ∈ p.returnvals)
            push!(returnvals,figh)
        end
        #if(p.savedata)
        #    SavePlots(figh,p.path,"Tmap")
        #end
        if(p.savedata)
            mkdelim(p.path * "transmission.txt", [p.E_samples TofE])
            mkdelim(p.path * "Tmap.txt", [Tmap'])
            #mkdelim(p.path * "Tmap.txt", [vec(S.*kgrid) vec(Tmap')])
            #jldsave(p.path * "fig_Tmap.jld2", figh)
            #jldsave(p.path * "fig_transmission.jld2", figT)
            #SavePlots(figT,p.path,"transmission")
        end
        return TofE	
end

# takes in the parameters defined in runs.jl, returns H(k)
function runSCF(p,A,returnvals,save=false,path="",name="")
	H = Hgen(p,A,returnvals)
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
        returnvals = [] # array of things to return
        if(p.savedata)
            mkpath(p.path)
            mktxt(p.path * "params.txt",string(p))
        end
        println("=========== MWSM devices ==============")
        println("performing: $(p.sweep)")
        if(p.verbose)
            println("Running simulation with parameters: ")
            println("params = $p") # oh this is so ugly, clean it later
            if(save) println("Saving outputs to $path\n") end

            println("\n")
        end
        if(p.bands)
            println("Generating static hamiltonian...")
            H = runSCF(p,A,returnvals,save,path)
            testH = H([0;0;0])
            println("size of H = $(size(testH))")
            Q = I(p.n)⊗I(p.nsite)⊗τ₁⊗σ₀
            #Q = distFromDW(p,RvalsGen(p))⊗I(p.norb)⊗(2) 
            #Q = zpos(RvalsGen(p))⊗I(p.norb)⊗I(2)
            runBands(p,2^7,H,Q,true,p.arpack)
        end
        #println("testH₀ = \n $testH")
	# some arbitrary operators for bands routine
	
	# project onto Q = |In><In| for bands purposes
	# [unit cell] ⊗ [A/B site] ⊗ [atom type] ⊗ [px, py] ⊗ [spin]
	#DOS, Evals = runDOS(20,H,λ,save,path,Beff)
	#runLDOS(20, H, λ,save,path,true,Beff)
        if(p.transport)
            TofE = NEGF_2contacts_1layer(p,A,returnvals)
            println("Done!\n")
            #push!(returnvals,TofE)
            #return returnvals
            return TofE
        else
            return returnvals
            println("Done!\n")
        end
end
end
