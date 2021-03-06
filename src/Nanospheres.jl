

#push!(LOAD_PATH, "./src/")
module Nanospheres

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

export peaksNS, NSprofileGen, pseudoBstats, closestPeak, strain, testhess, hoppingModification

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


function hoppingModification(pureNNs, λ, dux_dx, dux_dy, duy_dx, duy_dy, dh_dx, dh_dy, Ax, Ay)
	NNs = deepcopy(pureNNs)
	β = 3.37

	# loop over all of the nearest-neighbor bonds
	for NN in NNs
		gA = NN.gA # grid-value of site 1

		# approximate bond as halfway between NNs, from A site
		R = radius([1;1],gA,λ) + (1/2)*NN.r
		dx = NN.r[1]; dy = NN.r[2];

		# generate x,y,z components of modified bonds 
		Δh = dh_dx(R)*dx + dh_dy(R)*dy;
		Δux = dux_dx(R)*dx + dux_dy(R)*dy;
		Δuy = duy_dx(R)*dx + duy_dy(R)*dy;
		

		#total length of strained bond
		L = √(Δh^2 + (dx+Δux)^2 + (dy+Δuy)^2)
		#println("dx, dy, dz = $Δh, $Δux, $Δuy. Length = $L")
		# we will apply the relevant phase shift from the gauge field, exp(i*2π/Φ₀*∫A⋅δL)
		Φ₀ = h/q;
		ϕ = 2*π/Φ₀
		
		# integrate -- approximate ∫A⋅dL = 𝐀⋅𝛅
		AdL = Ax(R)*(dx+Δux) + Ay(R)*(dy+Δuy)
		ϕₚₛ = ϕ*AdL	

		#apply the phase shift: 
		#NN.t = t*exp(im*ϕₚₛ)
		NN.t = t*exp(-β*(L/δ - 1))*exp(im*ϕₚₛ)
		#println("Hopping integral t = $(NN.t)")

	end
	return NNs
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
