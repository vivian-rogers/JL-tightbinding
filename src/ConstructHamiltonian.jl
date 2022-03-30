

module ConstructHamiltonian
#push!(LOAD_PATH, "./src/")
using LinearAlgebra
using SparseArrays
using Operators
using Arpack
using MaterialParameters
using UsefulFunctions
using Printf
using Constants
using nearestNeighbors
#using PlotStuff
using GenNNs

export Hgen, makeH


function Hgen(spin)

	# modify the hopping and generate the static part of the hamiltonian
	#NN2s = nn2(λ,Rvals,U₂)
	# debug
	#show(H_hop)
	#println("edgeNNs")
	#show(edgeNNs)

	#Himp = chargeImpurities(λ,fz,R,V₀)
	
	# should be imported, is estimate for now
	H_onsite = Diagonal(ϵ_onsite)
	
	NNAtoms = NNatoms()
	NN_hops = NNfromAtoms(NNAtoms)
	H_hop, edge_hops = nnHoppingMat(NN_hops,4)
	H₀ = H_onsite .+ H_hop # add on Hpot, Hcharge, H+U, etc here
	#H₀ = sparse(H_hop .- μ*I(N) .+ H_U₂ .+ Himp) # add on Hpot, Hcharge, H+U, etc here
	
	#=function H_pre(k,B) #pre-scf
		H_edge = spzeros(ComplexF64, N,N)
		for NN in edgeNNs
			H_edge[NN.a,NN.b] += NN.t*exp(im*k⋅(NN.N[1]*a₁+NN.N[2]*a₂)) 
		end
		return dropzeros(H_edge.+Hstatic)
	end=#
	
	# this may be bad, take out later
	#B = 5
	#H_U = Diagonal(deepcopy(U*(min.(n₊,n₋)))) # set the avg energy around 0
	#H_B = Diagonal(deepcopy(B*(n₊.-n₋)))
	#function H(k::Array{Int64,1})

	N = size(H_onsite)[1]
	
	β = 0.6*eV
	function H(k)
		H_edge = zeros(ComplexF64, N,N)
		
		for NN in edge_hops
			H_edge[NN.a,NN.b] += NN.t*exp(im*k⋅(A*NN.N)) 
		end
		if(spin == true)
			return Hermitian( (H_edge .+ H₀)⊗I(2) .+ β*I(2)⊗[1 0; 0 0]⊗σ₃)
		else
			return Hermitian(H_edge .+ H₀)
		end
		#return I(4)
	end
	println("Great, SCF converged. Returning H(k).\n")
	return H 
	#return H, n₊, n₋, R, pseudoB.(R)
	#return H, n₊, n₋, R, pseudoB.(R)
end


#function makeH(H₀, edgeNNs, λ, k::Array{Int64,1})
function makeH(H₀, edgeNNs, λ, k)
	
	# a jank thing to make the compiler happy, speeds up runtime × λ² = 8000 for real system
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
