

#push!(LOAD_PATH, "./src/")
module VectorPotential

#using Plots
#using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using SparseArrays
using Operators
using GMaterialParameters
using Bands
using ForwardDiff
using nearestNeighbors
#using Interpolations

export peaksNS, NSprofileGen, pseudoBstats, closestPeak, Bvals, strain, testhess, hoppingModification, checkPeriodicField, zeeman

function closestPeak(λ)
	Rarray = Rvals(λ) 
	RvalsNS = peaksNS(λ)
	N = 2*λ^2
	peakDist = zeros(N)
	for i in 1:N
		r1 = Rarray[i,:]
		#show(RvalsNS)
		Δx = minimum(map(r2->norm(r1-r2),RvalsNS))
		peakDist[i] = Δx[1]
	end
	return sparse(Diagonal(peakDist))
end

function pseudoBstats(Bvals)
	Bmax = maximum(Bvals)
	Bavg = avg(abs.(Bvals))
	Brms = √(avg(Bvals.^2))
	return Brms, Bavg, Bmax 
end


function peaksNS(λ,θ=0) #generates the positions of the nanospheres at each of the corners of the unit cell
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	#needs work to rotate
	#RvalsNS = Array[]
	#RvalsNS = Array{Array{}(1),}(1)
	#RvalsNS = Array{Array{Float64, 2}}(4)
	#RvalsNS = Array{::Array{Float64,1},1}(undef, 4)
	RvalsNS = zeros(4,2)
	for i = 0:1; 
		for j = 0:1
			x = i*a₁ + j*a₂
			RvalsNS[2*j+i+1,:] = x
			#push!(RvalsNS,x)
		end
	end
	return [RvalsNS[i, :] for i in 1:size(RvalsNS, 1)]
end

#function height(Δh0,σ,RvalsNS,x,y)
#	return sum([Δh0*Gaussian(σ,[x,y],RvalsNS[i,:]) for i = 1:4]);

function rotate(θ)
	return [cosd(θ) -sind(θ); sind(θ) cosd(θ)]
end

function NSprofileGen(λ,Δh0,σ,θ=0) # generates the z(x,y) profile of the nanospheres
	Δedge = 2*nm
	#Δx = 0.05*nm
	nGrid = 10*λ
	RvalsNS = peaksNS(λ)
	MaxR = maximum(RvalsNS[4][:])

	# evenly spaced grid -- not used for much
	xVals = range(0 - Δedge, stop=MaxR + Δedge, length=nGrid);
	yVals = range(0 - Δedge, stop=MaxR + Δedge, length=nGrid);
	
	# this function important -- passed around everywhere
	a₁ = secd(θ)*λ*a*[1;0]; a₂ = secd(θ)*λ*a*[cosd(60);sind(60)]; 
	#z(R) = Δh0*Gaussian(σ,R,((1/2)*a₁ + (1/2)*a₂))
	z(R) = sum([Δh0*Gaussian(σ,R,rotate(θ)*(i*a₁ + j*a₂)) for i = -2:3 for j = -2:3]);
	#z(R) = Δh0*sin(2*π*R⋅(a₁+a₂).^(-1));
	zVals = [z([x,y]) for y in yVals, x in xVals];
	#surface(xVals,yVals,z)
	return xVals, yVals, zVals, z
end

function Gaussian(σ,R,Rcenter) # generate a gaussian around Rcenter
	C = 1
	#show(Rcenter)
	return Float64(C*exp((-1/2)*norm(R - Rcenter)^2/σ^2))
end	

function testhess(z)
	dux_dx(R) = a*Hessian(z,R)[1,1]
	return dux_dx
end


function hoppingModification(pureNNs, A::Function)
	NNs = deepcopy(pureNNs)
	# loop over all of the nearest-neighbor bonds
	AofR = zeros(3)
	for NN in NNs
		# approximate bond as halfway between NNs, from A site
		R = NN.ra + (1/2)*NN.r
		Φ₀ = h/q;
		ϕ = 2*π/Φ₀
		
		# integrate -- approximate ∫A⋅dL = 𝐀⋅𝛅
		AofR = A(R)
		#apply the phase shift: 
		∫AdL = AofR⋅NN.r
		#NN.t = t*exp(im*ϕₚₛ)
		NN.t = exp(im*ϕ*∫AdL)*NN.t
	end
	return NNs
end

function B(A::Function,R::Vector{Float64})
	return curl(A,R)
end

function checkPeriodicField(A::Function, p, npts::Int=30)
	for i = 1:npts
		δ = rand()*p.SLa₁ + rand()*p.SLa₂ + rand()*p.SLa₃
		B₀ = B(A,δ)
		cutoff = sum(B(A,δ).^2)*10^-4
		n = 1
		for ix = -n:n
			for iy = -n:n
				B₁ = B(A,δ .+ ix*p.SLa₁ .+ iy*p.SLa₂)
				if(sum((B₀ - B₁).^2) > cutoff)
					println("B₀ = $(round.(B₀,sigdigits=4)), B₁ = $(round.(B₁,sigdigits=4))")
					println("B Field not periodic! Check A(R) definition.")
					println("Disregard if modelling finite device")
					#return false
				end
			end
		end
	end
	println("Gauge potential appears periodic with superlattice.")
	println("Net flux may still be nonzero, should consider checking this.")
	println("May proceed with A(R) for superlattice modelling.")
	return true
end

function curl(f::Function, R::Vector{Float64},δ::Float64 = 10^-14)
	δ1 = [1;0;0]*δ; δ2 = [0;1;0]*δ; δ3 = [0;0;1]*δ;
	curlx = f(R+δ2)[3] - f(R-δ2)[3] - f(R+δ3)[2] + f(R-δ3)[2]
	curly = f(R+δ3)[1] - f(R-δ3)[1] - f(R+δ1)[3] + f(R-δ1)[3]
	curlz = f(R+δ1)[2] - f(R-δ1)[2] - f(R+δ2)[1] + f(R-δ2)[1]
	return [curlx; curly; curlz]*(2*δ)^-1
end

function Bvals(A::Function, Rvals::Vector{Vector{Float64}})
	return curl.(A,Rvals)
	#return B = [curl(A,R) for R in Rvals]
end

function zeeman(Bvals::Vector{Vector{Float64}},  p)
	# only defined for S-like orbitals with lz = 0
	N = p.n*p.nsite*p.norb*2
	#for i = 1:N
	zeeman = spzeros(ComplexF64, N, N)
	C = ħ/(2*m₀) #sans q factor -> eV
	for i in eachindex(Bvals)
		site = zeros(p.n*p.nsite*p.norb); site[i] = 1
		B = Bvals[i]
		#show(size(zeeman))
		#println("now size of additional site")
		#show(size(site⊗σ₁))
		zeeman .+= 2*C*Diagonal(site)⊗(B[1]*σ₁ .+ B[2]*σ₂ .+ B[3]*σ₃)
	end
	return zeeman
end

function strain(λ,Δh0,fᵤ,σ,z)
	
	# this function full of function definitions -- start defining the functions 
	
	∇h(R) = ∇(z,R)
	dh_dx(R) = ∇h(R)[1]
	dh_dy(R) = ∇h(R)[2]


	# generate in-plane displacement profile
	u = R -> -fᵤ*a*∇(z,R)
	# we don't actually need this though, just the derivatives
	# little bit jank but it gets the job done 
	dux_dx(R) = -fᵤ*a*Hessian(z,R)[1,1]
	duy_dx(R) = -fᵤ*a*Hessian(z,R)[1,2]
	dux_dy(R) = -fᵤ*a*Hessian(z,R)[2,1]
	duy_dy(R) = -fᵤ*a*Hessian(z,R)[2,2]

	# as defined in ENGINEERING RELATIVISTIC FERMIONS IN CONDENSED MATTER SYSTEMS, Youngseok Kim thesis 
	ε_xx = R -> dux_dx(R) + (1/2)*(dh_dx(R))^2
	ε_yy = R -> duy_dy(R) + (1/2)*(dh_dy(R))^2
	ε_xy = R -> (1/2)*(dux_dy(R)+duy_dx(R)) + (1/2)*(dh_dx(R)*dh_dy(R))
	ε_eff = R -> ε_xx(R) + ε_yy(R)
	κ = 1; β = -3.37; Φ₀ = h/q;
	C = κ*β*Φ₀/(a^2*4*π)
	Ax = R -> C*a*(ε_xx(R)-ε_yy(R))
	Ay = R -> -C*a*(2*ε_xy(R))
	pseudoB(R) = ∇(Ay,R)[1] - ∇(Ax,R)[2] # this is the curl

	return pseudoB, ε_eff, dux_dx, dux_dy, duy_dx, duy_dy, dh_dx, dh_dy, Ax, Ay
end

end
