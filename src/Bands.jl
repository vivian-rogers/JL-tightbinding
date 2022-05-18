

module Bands
using UsefulFunctions
#using Arpack
using LinearAlgebra
using Arpack
#using Suppressor
using Logging
#IJulia.installkernel("Julia nodeps", "--depwarn=no")

export getBands, project, expectedValue

function interpolate(n,klist,a,kdict)
	#=kdict = Dict(
		    "Γ" => [0;    0;    0],
		    "X" => [1/2;  0;    0],
		    "K" => [√2/4; √2/4; 0],
		    "U" => [-1/2; 0;   0],
		    "A" => [0;   0; 1/2],	
		    "H" => [1/3;1/3;1/2]
		    )
	=#
	#=kdict = Dict(
		    "Γ" => [0;   0;   0],
		    "K" => [1/2; tand(30)/2; 0],
		    "M" => [1/2; 0;   0],
		    "M-" => [-1/2; 0;   0],
		    "A" => [0;   0; 1/2],	
		    "H" => [1/3;1/3;1/2]
		    )=#	    
	#highSym = map(k->kdict[k], klist)
	highSym = (1)*map(k->kdict[k], klist)
	kpts = Any[]
	d = 1/n
	for k1 in 1:(size(highSym)[1]-1)
		k2 = k1 + 1
		for f = 0:d:(1-d) #fraction of the 2nd kpt
			k = (1-f)*highSym[k1] + f*highSym[k2]
			push!(kpts,k)
		end
	end
	push!(kpts,highSym[end])
	return kpts
end

function getBands(klist, kdict, n, a, Hofk) #takes in array of kpt strings, number of interpolation pts, H(k) function
	kpts = interpolate(n,klist,a,kdict)
	nk = size(kpts)[1]
	
	#	show(kpts)
	testH = Hofk([0;0;0])
	if(any(isnan,testH)==true)
		throw(DomainError(testH, "Something broken in hamiltonian definition! Returning NaN"))
		return
	end
	#initialize the band array
	#λ_test, evecs_test = eig(testH, 1E-12)
	arpack = false # small hamiltonian, few bands
	if(arpack)
		maxiter = 400
		nE = 4*Int(floor(log2(size(testH)[1])))
		nEig = size(testH)[1]
		if(nE < size(testH)[1])
			println("Heads up! $nE / $nEig eigvls are being calculated")
		end
		λ_test, evecs_test = eigs(testH,nev=nE, maxiter=maxiter)
		#nE = size(λ_test)[1]
		ndim = size(evecs_test)[2]
	else
		nE = size(testH)[1]
		nEig = nE
	end
	Evals = zeros(Float64, nk, nE)
	Estates = zeros(ComplexF64, nk, nE, nE) 
	#go through and fill it up
	
	#
	#d = 100
	#Logging.disable_logging(Logging.Warn)
	for ik in 1:nk
		k = kpts[ik]
		H = Hofk(k)
		#print("\n\nnk = ")
		#show(k)
		#print("\nH(k) = ")
		#show(H)
		#Eofk, Estatek = eigs(Hermitian(H))
		#Eofk, Estatek = eigen(H)
		if(arpack)
			Eofk, Estatek = eigs(H,nev=nE, which=:SM, maxiter=maxiter)
		else
			Eofk, Estatek = eigen(H)
		end
		#show(size(Estatek))
		#Eofk, Estatek = eigen(H)
		#Eofk = eigvals(H)
		for iE in 1:nE
			Evals[ik,iE] = real(Eofk[iE]) 
		end
		#Estatesk = eigvecs(H)
		for iE1 in 1:nE
		#for iE1 in 1:nEig; for iE2 in 1:nEig
			Estates[ik,:,iE1] = Estatek[:,iE1] # pre-arpack
				#Estates[ik,iE1,iE2] = Estatek[iE2,iE2]
		end
	end
	return Evals, Estates
end

function project(projKet, Estates) #takes in array of kpt strings, number of interpolation pts, H(k) function
	nk = size(Estates)[1]
	nN = size(Estates)[2]
	nE = size(Estates)[3]
	projVals = zeros(Float64, nk, nE)
	adjoint = projKet' #takes |proj> -> <proj|
	for ik in 1:nk
		for iE in 1:nE
			projVals[ik,iE] = abs(adjoint*Estates[ik,iE,:])^2 #performs <proj|band>
		end
	end
	return projVals
end

function expectedValue(Q, Estates) #takes in array of kpt strings, number of interpolation pts, H(k) function
	nk = size(Estates)[1]
	nN = size(Estates)[2]
	nE = size(Estates)[3]
	projVals = zeros(Float64, nk, nE)
	#adjoint = projKet' #takes |proj> -> <proj|
	for ik in 1:nk
		for iE in 1:nE
			#show(size(Estates))
			Estate = Estates[ik,:,iE]
			projVals[ik,iE] = real(Estate'*Q*Estate) #performs <proj|band>
		end
	end
	return projVals
end

#function plotBands(klist, n, E)
#	nk = length(klist)
#	indices = LinRange(0,nk,nk*n+1) 
#	
#
#end

end
