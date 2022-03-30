

module DOS
using UsefulFunctions
#using Arpack
using LinearAlgebra
using Arpack
using GMaterialParameters
using Logging
#IJulia.installkernel("Julia nodeps", "--depwarn=no")

export getDOS, getLDOS

function getLDOS(kgrid, weights, Hofk, λ, Emax=4, Emin=-4, Ebins=100) #takes in array of kpt strings, number of interpolation pts, H(k) function
	testH = Hofk([0;0])
	nk = size(kgrid)[1]
	#initialize the band array
	#λ_test, evecs_test = eig(testH, 1E-12)
	N = 2*λ^2
	maxiter = 400
	nE = 12*Int(floor(log2(size(testH)[1])))
	nEig = size(testH)[1]
	if(nE < size(testH)[1])
		println("Heads up! $nE / $nEig eigvls are being calculated")
	end
	λ_test, evecs_test = eigs(testH,nev=nE, maxiter=maxiter)
	ndim = size(evecs_test)[2]
	LDOS = zeros(Ebins,2*λ+1)
	DOS = zeros(Ebins)
	Evals = LinRange(Emin,Emax,Ebins)
	function idos(E)
		if(E < Emax && E>Emin)
			return Int(round((E-Emin)*Ebins/(Emax-Emin)))
		else
			return false
		end
	end
	Logging.disable_logging(Logging.Warn)
	for ik in 1:nk
		k = kgrid[ik]
		w = weights[ik]
		H = Hofk(k)
		Eofk, Estatek = eigs(H,nev=nE, which=:SM, sigma=10^-5, maxiter=maxiter)
		for iE in 1:nE
			Eval = real(Eofk[iE]) # particular energy value
			Estate = Estatek[:,iE]
			index = idos(Eval)
			if(index != false)
				# add to real DOS
				DOS[index] += w
				#take a diagonal path!
				for i = 1:λ
					i1 = g2i(i,i,λ)
					overlap = abs(Estate[i1])^2
					LDOS[index,2*i-1] += w*overlap
				end
				for i = 1:λ
					i2 = g2i(i+1,i,λ)
					overlap = abs(Estate[i2])^2
					LDOS[index,2*i] += w*overlap
				end
			end
		end
	end
	return LDOS, DOS, Evals
end
function getDOS(kgrid, weights, Hofk, Emax=1, Emin=-1, Ebins=100) #takes in array of kpt strings, number of interpolation pts, H(k) function
	testH = Hofk([0;0])
	nk = size(kgrid)[1]
	#initialize the band array
	#λ_test, evecs_test = eig(testH, 1E-12)
	maxiter = 400
	nE = 4*Int(floor(log2(size(testH)[1])))
	nEig = size(testH)[1]
	if(nE < size(testH)[1])
		println("Heads up! $nE / $nEig eigvls are being calculated")
	end
	λ_test, evecs_test = eigs(testH,nev=nE, maxiter=maxiter)
	ndim = size(evecs_test)[2]
	DOS = zeros(Ebins)
	Evals = LinRange(Emin,Emax,Ebins)
	function idos(E)
		if(E <= Emax && E>=Emin)
			return Int(round((E-Emin)*Ebins/(Emax-Emin)))
		else
			return false
		end
	end
	
	# I do this so that I am not yelled at by the eigensolver for shifting the spectrum
	Logging.disable_logging(Logging.Warn)
	for ik in 1:nk
		k = kgrid[ik]
		w = weights[ik]
		H = Hofk(k)
		Eofk, Estatek = eigs(H,nev=nE, which=:SM, sigma=10^-5, maxiter=maxiter)
		for iE in 1:nE
			Eval = real(Eofk[iE]) # particular energy value
			index = idos(Eval)
			if(index != false)
				DOS[index] += w
			end
		end
	end
	return DOS, Evals
end

end
