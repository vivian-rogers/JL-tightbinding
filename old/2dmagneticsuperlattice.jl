push!(LOAD_PATH, "./src/")

using Plots
using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
#using GrapheneParameters
using Bands



⊗(A,B) = kron(A,B)

function rotate(θ::Float64)
	return [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
end

a = 1
function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	m = 0*eV # mass term at M point
	t₁ = 1.3*eV # hopping parameter
	t₂ = 0.4*eV
	t₃ = 0.2*eV
	α = 0.01*eV
	α₂ = 0.0*eV
	ε = 0
	H₀ = zeros(ComplexF64,2,2)
	Hᵣ = zeros(ComplexF64,4,4)
	A = [0; 0.4; 0] #not actually the gauge field but close enough
	β = [0; 0; A[2]/2]*eV
	# 2nd nearest neighbors
	Rᵢ = a*[[1;0;0],[-1;0;0],[0;1;0],[0,-1,0]]
	for R in Rᵢ
		#δᵢ = rotate((π/2)*i)*δ
		H = -t₂*exp(im*(k⋅R + A⋅R))*I(2)
		H₀ = H₀ .+ H
		soc = (R[1]*σ₂ .- R[2]*σ₁)
		Hᵣ .+= α₂*exp(im*(k⋅R))*I(2)⊗soc
	end
	# first nearest neighbors
	aᵢ = a*[[1;0;0],[1;1;0],[0;1;0]]
	δᵢ = (a/2)*[[1;-1;0],[-1;1;0],[1;1;0]]
	
	# over the xy edges
	δ = (a/2)*[-1;-1;0]
	H₀ .+= t₁*[    0    exp(im*(δ⋅A))
		    exp(-im*(δ⋅A)) 0]
	for i in eachindex(δᵢ)
		R = aᵢ[i]
		δ = δᵢ[i]
		#δᵢ = rotate((π/2)*i)*δ
		H = t₁*[	
		       0 		exp(im*(R⋅k + δ⋅A));
		       exp(-im*(R⋅k + δ⋅A))  0
		       ]
		H₀ = H₀ .+ H
	end
	
	# construct SoC term
	
	aᵢ = a*[[1;0;0],[1;1;0],[0;1;0]]
	δᵢ = (a/2)*[[1;-1;0],[-1;1;0],[1;1;0]]
	
	# in-unit cell
	δ = (a/2)*[-1;-1;0]
	#Hᵣ .+= α*-im*τ₂⊗(δ[1]*σ₂ - δ[2]*σ₁)
	Hᵣ .+= α*τ₁⊗(δ[1]*σ₂ .- δ[2]*σ₁)
	# over the xy edges
	for i in eachindex(δᵢ)
		δ = δᵢ[i]
		R = aᵢ[i]
		
		# add SoC term
		soc = (δ[1]*σ₂ .- δ[2]*σ₁)
		Hᵣ .+= α*[0 exp(im*k⋅R)
			  exp(-im*k⋅R) 0]⊗soc
	end
	# over z
	#=aᵢ = a*[[0;0;1],[0;0;-1]]
	δᵢ = aᵢ
	for i in eachindex(δᵢ)
		δ = δᵢ[i]
		R = aᵢ[i]
		H₀ .+= t₃*exp(im*k⋅R)*I(2)
		# add SoC term
		soc = (δ[3]*σ₁ - δ[3]*σ₂)
		Hᵣ .+= exp(im*k⋅R)*I(2)⊗soc
	end=#
	H = H₀⊗I(2) .+ Hᵣ .+ β[3]*I(2)⊗σ₃
	return Hermitian(H)
end

function B(R)
	B₀ = [0,0,5]
	return B₀
end



#=function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	#Hᵣ = α*( km[1]*km[2]^2*σ₁ - km[2]*km[1]^2*σ₂ )
	#Hᵣ = α*( sin(k[1])*sin(k[2])^2*σ₁ - sin(k[2])*sin(k[1])^2*σ₂ )
	#Hᵣ = α*(σ₁*km[2] .- σ₂*km[1])*√(abs(1-km[1]-km[2]))
	#Hᵣ = α*(σ₁*km[2] - σ₂*km[1])
	#Hᵣ = α*I(2)⊗(σ₂*km[3] - σ₃*km[2] + σ₃*km[1] - σ₁*km[3] + σ₁*km[2] - σ₂*km[1])

	β = [0,0,0.01]*eV # magnetization vector
	α = 0.3*eV # SoC term
	m = 0*eV # mass term at M point
	t = 1.0*eV # hopping parameter
	
	km = [cos(k[1]), cos(k[2]), sin(k[3])]
	# basic model of square lattice bands
	H₀ = 0*( -t*(km⋅km)*I(2) .+ 3*t*km[1]*km[2]*(1-km[3])*τ₃ )⊗I(2)
	
	# add τ₁⊗(𝐤⋅𝛔) dirac cone term
	Hᵣ = α*τ₁⊗(km[1]*σ₁ .+ km[2]*σ₂ .+ km[3]⊗σ₃)
	#Hᵣ = α*τ₁⊗(km[1]*σ₁ .+ km[2]*σ₂ .+ km[3]⊗σ₃)*√(abs(1-km[1]-km[2]))
	
	# Mass term
	#Hₘ = m*τ₃⊗I(2) 
	
	# Zeeman splitting + SoC cross term
	#Hᵦ = (β[1]⊗σ₁ .+ β[2]⊗σ₂ .+ β[3]⊗σ₃)
	Hᵦ = β[3]*I(2)⊗σ₃ .+ τ₃⊗(β[1]*σ₁.+β[2]*σ₂)
	
	H₀ = zeros(ComplexF64,2,2)
	for i = 1:3
		for d = [-1,1]
			r = zeros(3)
			r[i] = d/2
			#Hi = zeros(
			H₀ .+= d*exp(im*k⋅r)*(r[1]*k[1]*σ₁ .+ r[2]*k[2]*σ₂ .+ r[3]*k[3]*σ₃)
			#H₀ .+= exp(im*k⋅r)*(r[1]*σ₁ .+ r[2]*σ₂ .+ r[3]*σ₃)
			#H[2,1] += exp(-im*k⋅r)
			#H .+= τ₁
		end
	end
	H  = τ₁⊗H₀ .+ Hᵦ
	return H
end =#


klist = ["Γ","M","X₁","M","X₃","M","X₂","Γ"]
klist = ["Γ","M","X₁","Γ"]
C = 2*π
kdict = Dict(
	    "Γ" => C*[0;   0;   0],
	    "R" => C*[1/2; 1/2; 1/2],
	    "X₁" => C*[1/2; 0;   0],
	    "X₂" => C*[0; 1/2;   0],
	    "X₃" => C*[0; 0;   1/2],
	    "M" => C*[1/2; 1/2;   0],
	    )	    
nk = 2048
println("Getting eigenvalues of 2D weyl lattice between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, kdict, nk, a, H)
#display(27.2*E)
println("Plotting...")
plotBands(klist,nk,E)
#plot3DSpectra(H,kdict["M"]+[0;0;0.00*π/a],a,160,160,0.1,0.1)
#plot3DSpectra(H,kdict["M"],a,160,160,0.1,0.1)
println("Done! Press ctrl+d to quit")

