module TB_GD_Helper
using LinearAlgebra
using SparseArrays

export argMinimize

#function argMinimize(f::Function, startingR, )
	




function tbErrorGen(EvsK_DFT, kpts, nE = 5)
	# read in bands file from DFT
	function tbError(H)
		err = 0
		nk = size(kpts)[1]
		for ik in 1:nk
			fullEofk, fullEigstatek = eigen(H(k))
			Eofk_DFT = [fullEofk[i] = 1:nE]
			Eofk_TB = EvsK_DFT[i]
			for iE = 1:nE
				err += (EofK_DFT[iE] - Eofk_TB[iE])^2
			end
		end
		return √(err) / (nE*nk)
	return tbError
end

function load_DFT_bands(filename)
	# 
	readdlm(filename, '\t', Int, '\n')
	
end

function totalIndex(params,atom,orbital,spin=0)
	iA = params.atomicEigenBasis[atom] - 1
	iO = params.orbitalEigenBasis[orbital] - 1
	iS = 0
	if(spin == 0):
		# return itot, iA, iO
		return (iA*iO + iO + 1), iA+1, iO+1
	else:
		if(spin == "up"):
			iS = 1
		elseif(spin == "down"):
			iS = 0
		end
		# return itot, iA, iO, iS
		return (2*(iA*iO + iO + 1) + iS), iA+1, iO+1, iS+1
	end
end

function GD_Hinit_ScN(param, coeffs)
	# generate list of 
	# Return H₀, Coeffs, Xindices, Yindices, phaseOffset in C*exp(val*k)
	nA = size(param.atoms)[1] # size of atom basis
	nO = 9 # size of orbital eigenbasis
	nspin = 2 # account for spin (1: no spin at all, 2: hopping independent of spin, 3: full SoC)
	n = nA*nO
	H₀ = zeros(n,n)
	i = 1 # index of coeff
	coeffs = []; Xinds = []; Yinds = []; POs = [];
	# first coeffs determine onsite energies
	for atom in param.atoms:
	        for orbital in atom.orbitals:
			iH, iA, iO = totalIndex(param,atom,orbital)
			#orbitalOn = atom.orbitals[iO]
			#if(orbitalOn == 0):
			H₀[iH,iH] = coeffs[i]
			i += 1
			#end
		end
	end
	
	H₀ = H₀⊗I(2) # put into spin eigenbasis

	# then determine hopping interactions
	
	for bondtype in param.connections:
		Atom1 = bondtype.Atom1;
		Atom2 = bondtype.Atom2;
		scale = param.radiusCoeff
		NNbonds = scale.*param.NNs[bondtype.NN]
		for bond in NNbonds:
			for orbital1 in Atom1.orbitals:
				iH1, iA1, iO1 = totalIndex(param,Atom1,orbital1)
				for orbital2 in Atom2.orbitals:
					iH2, iA2, iO2 = totalIndex(param,Atom2,orbital2)
					t = coeffs[i]
					
					i += 1
	
end

function GD_Hbuilder(H₀, Coeffs, Iinds, Jinds, phaseOffsets)
	function H(k) 
		return H₀ .+ sparse(Iinds,Jinds,[PO⋅k for PO in phaseOffsets])
	end
	return H
end

function TB_GD_Wrapper(nk, klist, kdict, params, DFT_filename::String)

	DFT_bands = load_DFT_bands(DFT_filename)
	H₀, coeffs, Iinds, Jinds, POs = GD_Hinit_ScN
	H = GD_Hbuilder(H₀, coeffs, Iinds, Jinds, POs)
	

	#	show(kpts)
	testH = Hofk([0;0;0])
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
	for i in eachindex(kpts)
		k = eachindex[i]
		EH(k)
	  



