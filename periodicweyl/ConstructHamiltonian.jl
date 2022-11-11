
push!(LOAD_PATH, "./src/")
push!(LOAD_PATH, "../src/")

module ConstructHamiltonian
#push!(LOAD_PATH, "./src/")
using LinearAlgebra
using SparseArrays
using Operators
using Arpack
using MaterialParameters
using PlotStuff
using UsefulFunctions
using Printf
using Constants
using InBiNNs
using VectorPotential
using Distributions
using MakiePlots
#using PlotStuff
#using GenNNs

export Hgen, makeH

# Takes in the parameter list and the vector potential
function Hgen(p,A::Function,returnvals)

	# modify the hopping and generate the static part of the hamiltonian
	#NN2s = nn2(λ,Rvals,U₂)
	# debug
	#show(H_hop)
	#println("edgeNNs")
	#Himp = chargeImpurities(λ,fz,R,V₀)
	
	# should be imported, is estimate for now
	#H_onsite = Diagonal(ϵ_onsite)
	
	#NNAtoms = NNatoms()
	#NN_hops = NNfromAtoms(NNAtoms)
	NNs = genNNs(p)
	#Rvals = RvalsGen(p)⊗I(p.norb) 
	#println("NNs = $NNs")
	NNs = pruneHoppings(NNs,p.prune) # cut off the periodic boundary hoppings
	if(p.fieldtype == "A")
		NNs = hoppingModification(NNs,A) # apply the peierls phase
	end
	H₀, edge_NNs = nnHoppingMat(NNs,p) 
        NormalDist = Normal(0,p.μ_disorder)
        H_onsite = Diagonal(rand(NormalDist,p.n*p.nsite))⊗I(p.norb*2) .+ p.μ*I(p.n*p.nsite*p.norb*2)
        #H_onsite = p.μ_disorder*Diagonal(rand(p.n*p.nsite).-0.5)⊗I(p.norb*2) .+ p.μ*I(p.n*p.nsite*p.norb*2)
	#H_onsite = 0*3*p.t*I(p.n)⊗I(p.nsite)⊗τ₃⊗I(2)
	Rvals = RvalsGen(p)
        println("Applying magnetization profile...")
	# for plotting field at surface
	Rsurf = Vector{Float64}[]
	for Rval in Rvals
            if(Rval[3] ≈ (p.nz-1)*p.a₃[3])
                push!(Rsurf,Rval)
            end
	end
        if(p.deviceMagnetization==true)
            (Bfield, Bsurf, avgB) = fieldUtils(p,A,Rsurf,Rvals,returnvals)
            Hᵦ = zeeman(map(B -> Float64.(B), Bfield),p)
        else
            Hᵦ = 0I(p.n*p.nsite*p.norb*2)
        end
        #Bfield = fieldUtils(p,A,Rsurf)
	#println("Generating field")
	#Bsurf = zeros(size(Rsurf))
	#println("Plotting field at surface...")
	#plotScatter(Rsurf,Bvals(A,Rsurf)[3])
	#plotPoints(Rsurf,Bvals(A,Rsurf)[3])
	#plotFunct(Rsurf,Bvals(A)[3])
	#println("Calculating B field in slab...")
	#H₀ = 
	#println("B field = $Bfield T")
        #display(Hᵦ)
        #println("Zeeman splitting hamiltonian = $Hᵦ")
	#show(Hᵦ)
	#H₀ .+= Diagonal(rand(Float64,p.n*p.nsite*p.norb*2).-0.5)*10^-3
	#H₀ = H_onsite .+ H_hop # add on Hpot, Hcharge, H+U, etc here
	#H₀ = sparse(H_hop .- μ*I(N) .+ H_U₂ .+ Himp) # add on Hpot, Hcharge, H+U, etc here
	
	#H₀ = sparse(H₀ .+ Hᵦ)
	H₀ = sparse(H₀ .+ H_onsite .+ Hᵦ)
	#β = 0.02*eV*[1,0,0]
	#Hᵦ = I(p.nsite)⊗I(p.norb)⊗I(p.n)⊗(β[1]*σ[1] .+ β[2]*σ[2] .+ β[3]*σ[3])
	function H(k)
		if(any(isnan,k)==true)
				throw(DomainError(k, "Something broken in k vector definition! Returning NaN"))
				return
		end
		# hamiltonian describing the edges
		#Hₑ = zeros(ComplexF64, 2*p.nsite*p.norb*p.n,2*p.nsite*p.norb*p.n)
                #build sparse as Hₑ = sparse(rows, cols, elements)
                rows = Int[]; cols = Int[]; elements = ComplexF64[];
                for NN in edge_NNs
			
                        Δϕ = exp(im*k⋅(p.A*NN.N))
			#Δϕ = exp(im*k⋅(p.SLa₁*NN.N[1] + p.SLa₂*NN.N[2] + p.SLa₃*NN.N[3]))
			#=Hₑ[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]*Δϕ
			Hₑ[2*NN.b  , 2*NN.a-1] += NN.t[2,1]*Δϕ
			Hₑ[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]*Δϕ
			Hₑ[2*NN.b  , 2*NN.a  ] += NN.t[2,2]*Δϕ=#
                        for i = 1:2
                            for j = 1:2
                                push!(rows,2*NN.b+i-2);
                                push!(cols,2*NN.a+j-2);
                                push!(elements,copy(NN.t[i,j]*Δϕ))
                            end
                        end
		end
                Hₑ = sparse(rows,cols,elements)
		Htot = H₀ .+ Hₑ
                return dropzeros(Htot)
                #return dropzeros(Htot)
		#return dropzeros(Hᵦ)
		#return dropzeros(sparse(H₀ .+ H_onsite .+ Hₑ .+ Hᵦ))
		#return Hermitian(H₀ .+ H_onsite .+ Hₑ .+ Hᵦ)
	end
	println("Great, SCF converged. Returning H(k).\n")
	return H 
end


# all these functions are obsolete, but could be useful... especially chargeImpurities
# we are not adding any onsite +U terms quite yet, so the SCF loops may be unnecessary. 

function fieldUtils(p, A::Function, Rsurf::Vector{Vector{Float64}}, Rvals::Vector{Vector{Float64}}, returnvals)
	println("===== Calculating applied zeeman field properties =====")
	checkPeriodicField(A,p) # tells you if gauge field periodic with SL
	if(p.fieldtype == "A")
		println("Field type = vector potential, applying onsite zeeman and peierls term")
		Bfield = Bvals(A,Rvals); Nsites = p.n*p.nsite
		avgB = sum(Bfield)*Nsites^-1
		#avgB = [0;0;0]
		Bfield = [Bfield[i] .- avgB for i = 1:Nsites]
		
		Bsurf = Bvals(A,Rsurf)
		plot3Dvectors(Rvals,Bfield,[coeff*R[2] for R in Rvals],"x position (nm)", "y position (nm)", "β₂ (eV)")
		plotScatter(Rsurf,[B[1]-avgB[1] for B in Bsurf], "x position (nm)", "y position (nm)", "Bₓ (T)", "coolwarm",)
		println("ΣB field = $(sum(Bfield)) T")
		println("max B field = $(maximum(maximum.(Bfield))) T")
		return (Bfield, Bsurf, avgB)
		#return Bfield, Bsurf, avgB
		#println("Bsurf = $(Bsurf[:][:][3])")
		#Bfield = Bfield .- avgB
	elseif(p.fieldtype == "β")
		avgB = [0;0;0]
		coeff = ħ/m₀
		Bfield = A.(Rvals)
		Bsurf = A.(Rsurf)
                if(p.verbose)
                    println("Field type = exchange, applying onsite exchange-like term")
                    println("Σβ field = $(sum(coeff*Bfield)) eV")
                    println("max β field = $(coeff*maximum(maximum.(Bfield))) eV")
                end
                #plot3Dvectors(Rvals,Bfield,[coeff*B[2] for B in Bfield],"x position (nm)", "y position (nm)", "z position (nm)", "β₂ (eV)")
                if(p.plotfield)
                    colorfunc(B) = B⋅[0;1;0]
                    devicefig = render2TDevice(p,Rvals,coeff*Bfield,colorfunc,10.0) 
                    if("plotfield" ∈ p.returnvals)
                        push!(returnvals,devicefig)
                    end
                    #plotScatter(Rsurf,[coeff*B[2] for B in Bsurf], "x position (nm)", "y position (nm)", "β₂ (eV)", "coolwarm",)
                end
                #return (Float64.(Bfield), Float64.(Bsurf), avgB)
		return (Bfield, Bsurf, avgB)
	end
end
#function makeH(H₀, edgeNNs, λ, k::Array{Int64,1})
function makeH(H₀, edgeNNs, λ, k)
	
	#=let
		H_static = H_static_scf
		edgeNNs = edgeNNs_scf
	end =#

	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	N = 2*λ^2
	H_edge = spzeros(ComplexF64, N,N)
	for NN in edgeNNs
		#H_edge[NN.a,NN.b] += 1
		H_edge[NN.a,NN.b] += NN.t*exp(im*k⋅(NN.N[1]*a₁+NN.N[2]*a₂)) 
	end
	#show(sparse(H_edge .+ H_static))
	#return sparse(H_edge)
	return sparse(H_edge .+ H₀)
end


function chargeImpurities(λ, fz, R, V₀)
	N = 2*λ^2
	H = zeros(N)
	ζ = 1.5*nm
	zVals = fz.(R)
	#show(zVals)
	zs = zVals .- minimum(zVals)
	for i in 1:N # site to add voltage from
		R₁ = R[i]
		z = zs[i]*nm^-1 # add a factor of 10^9 if working in metres, this seems arbitrary 
		#show(z)
		for j in 1:N # site to add potential to
			R₂ = R[j]
			H[j] +=  z*exp(-norm(R₁-R₂)^2/(2*ζ^2))
		end
	end
	#show(H)
	return Diagonal(-V₀*H)
end

function testSCF(H, U, λ)
	println("========= Entering SCF loop ==========")
	N = 2*λ^2
	S = 2
	ne = S*λ^2 # number of electrons in system with μ = 0
	n₊ = rand(N) # confusing! lots of Ns, this n is the electron density vector
	n₋ = rand(N)# confusing! lots of Ns, this n is the electron density vector
	#n₋[1:2:N] .= 1
	eigsflag = false
	if(eigsflag)
		nev = Int(round(ne*1.1))
	else
		nev = N
	end
	err = 100; cutoff = 10^-7
	if(U == 0)
		println("no need for convergence, +U = 0")
		err = 0
	end
	#n[ne:N] = 1
	# generate an evenly-spaced k-grid here, maybe
	SCFd = 8
	println("SCF density: $SCFd")
	klist, kweights = kgrid(SCFd,λ)
	#kweights = [1] # weights on the k-grid
	nk = size(kweights)[1]
	loopn = 1
	B = 5
	while err > cutoff #SCF loop
		#print("SCF iteration $loopn: Error = $err")
		@printf "SCF iteration %i: Error = %g\n" loopn err
		
		# up spin calculation
		H_U = deepcopy(U*(min.(n₊,n₋))) # set the avg energy around 0
		n₊old = deepcopy(n₊)
		H_B = deepcopy(B*(n₊.-n₋))
		n₊ = zeros(N)
		for ik in 1:nk
			k = klist[ik]
			if(eigsflag)
				Htemp = sparse(H(k) .+ Diagonal(H_U) .+ Diagonal(H_B))
				Eofk, Estatek = eigs(Htemp, nev=nev, which=:SR)
			else
				Htemp = Hermitian(Array(H(k)) .+ Diagonal(H_U).+ Diagonal(H_B))
				Eofk, Estatek = eigen(Htemp)
			end
			occ_eigstate = fermi.(Eofk)
			for iE in 1:nev
				n₊ .+= real(kweights[ik]*occ_eigstate[iE]*(abs.(Estatek[:,iE])).^2)
			end
		end
		# down spin calculation
		n₋old = deepcopy(n₋)
		n₋ = zeros(N)
		for ik in 1:nk
		k = klist[ik]
			if(eigsflag)
				Htemp = sparse(H(k) .+ Diagonal(H_U) .- Diagonal(H_B))
				Eofk, Estatek = eigs(Htemp, nev=nev, which=:SR)
			else
				Htemp = Hermitian(Array(H(k)) .+ Diagonal(H_U).- Diagonal(H_B))
				Eofk, Estatek = eigen(Htemp)
			end
			occ_eigstate = fermi.(Eofk)
			for iE in 1:nev
				n₋ .+= real(kweights[ik]*occ_eigstate[iE]*(abs.(Estatek[:,iE])).^2)
			end
		end
		sum_n = sum(n₋ + n₊)
		@printf "Done: Δe⁻ / SL unit cell = %.2f - %.2f = %g\n" sum_n ne (sum_n -ne)
		nᵢ = ((sum_n-ne)/(λ^2*A_uc))*(cm/metre)^2
		P = avg(n₊) - avg(n₋)
		absP = avg(abs.(n₊ .- n₋))
		@printf "Total nᵢ = %g e⁻/cm²\n" nᵢ
		@printf "Total n₊ - n₋ = %g, abs(P) = %g\n\n" P absP
		err = norm(n₊-n₊old) + norm(n₋-n₋old)
		loopn += 1
	end
	return n₊, n₋
end
function SCF(H, U, λ)
	println("========= Entering SCF loop ==========")
	N = 2*λ^2
	S = 2
	ne = S*λ^2 # number of electrons in system with μ = 0
	n₊ = rand(N) # confusing! lots of Ns, this n is the electron density vector
	n₋ = rand(N)# confusing! lots of Ns, this n is the electron density vector
	#n₋[1:2:N] .= 1
	eigsflag = false
	if(eigsflag)
		nev = Int(round(ne*1.1))
	else
		nev = N
	end
	err = 100; cutoff = 10^-7
	if(U == 0)
		println("no need for convergence, +U = 0")
		err = 0
	end
	#n[ne:N] = 1
	# generate an evenly-spaced k-grid here, maybe
	SCFd = 0
	println("SCF density: $SCFd")
	klist, kweights = kgrid(3,λ)
	#kweights = [1] # weights on the k-grid
	nk = size(kweights)[1]
	loopn = 1
	while err > cutoff #SCF loop
		#print("SCF iteration $loopn: Error = $err")
		@printf "SCF iteration %i: Error = %g\n" loopn err
		
		# up spin calculation
		H_U₊ = U*(n₊.-1/2) # set the avg energy around 0
		n₊old = deepcopy(n₊)
		n₊ = zeros(N)
		for ik in 1:nk
			k = klist[ik]
			if(eigsflag)
				Htemp = sparse(H(k) .+ Diagonal(H_U₊))
				Eofk, Estatek = eigs(Htemp, nev=nev, which=:SR)
			else
				Htemp = Hermitian(Array(H(k)) .+ Diagonal(H_U₊))
				Eofk, Estatek = eigen(Htemp)
			end
			occ_eigstate = fermi.(Eofk)
			for iE in 1:nev
				n₊ .+= real(kweights[ik]*occ_eigstate[iE]*(abs.(Estatek[:,iE])).^2)
			end
		end
		#show(n₊)
		# down spin calculation
		H_U₋ = U*(n₋.-1/2) # set the avg energy around 0
		n₋old = deepcopy(n₋)
		n₋ = zeros(N)
		for ik in 1:nk
		k = klist[ik]
			if(eigsflag)
				Htemp = dropzeros(H(k) .+ Diagonal(H_U₋))
				Eofk, Estatek = eigs(Htemp, nev=nev, which=:SR)
			else
				Htemp = Hermitian(Array(H(k)) .+ Diagonal(H_U₋))
				Eofk, Estatek = eigen(Htemp)
			end
			occ_eigstate = fermi.(Eofk)
			for iE in 1:nev
				n₋ .+= real(kweights[ik]*occ_eigstate[iE]*(abs.(Estatek[:,iE])).^2)
			end
		end
		sum_n = sum(n₋ + n₊)
		@printf "Done: Δe⁻ / SL unit cell = %.2f - %.2f = %g\n" sum_n ne (sum_n -ne)
		nᵢ = ((sum_n-ne)/(λ^2*A_uc))*(cm/metre)^2
		P = avg(n₊) - avg(n₋)
		absP = avg(abs.(n₊ .- n₋))
		@printf "Total nᵢ = %g e⁻/cm²\n" nᵢ
		@printf "Total n₊ - n₋ = %g, abs(P) = %g\n\n" P absP
		err = norm(n₊-n₊old) + norm(n₋-n₋old)
		loopn += 1
	end
	return n₊, n₋
end



function fermi(E)
	T = 100.0 #temp in kelvin
	return (exp(E/(kB*T)) + 1)^-1
end

end
#function hoppingModification(
