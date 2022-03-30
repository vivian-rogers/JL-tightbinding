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
#= function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	m = 0*eV # mass term at M point
	t = 1.0*eV # hopping parameter
	H₀ = t*τ₁
	δ\_ = a*[0.5;0.5;0]
	for i = 1:3
		δᵢ = rotate((π/2)*i)*δ
		H = [	
		       0 		exp(im*(δᵢ⋅k));
		       exp(-im*(δᵢ⋅k))  0
		       ]
		H₀ = H₀ .+ H
	end
	return Hermitian(H₀)
end =#

function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	#Hᵣ = α*( km[1]*km[2]^2*σ₁ - km[2]*km[1]^2*σ₂ )
	#Hᵣ = α*( sin(k[1])*sin(k[2])^2*σ₁ - sin(k[2])*sin(k[1])^2*σ₂ )
	#Hᵣ = α*(σ₁*km[2] .- σ₂*km[1])*√(abs(1-km[1]-km[2]))
	#Hᵣ = α*(σ₁*km[2] - σ₂*km[1])
	#Hᵣ = α*I(2)⊗(σ₂*km[3] - σ₃*km[2] + σ₃*km[1] - σ₁*km[3] + σ₁*km[2] - σ₂*km[1])

	β = [0.005,-0.005,0.0]*eV # magnetization vector
	α = 0.3*eV # SoC term
	m = 0*eV # mass term at M point
	t = 1.0*eV # hopping parameter
	
	km = [cos(k[1]), cos(k[2]), sin(k[3])]
	# basic model of square lattice bands
	H₀ = ( -t*(km⋅km)*I(2) .+ 3*t*km[1]*km[2]*(1-km[3])*τ₃ )⊗I(2)
	
	# add τ₁⊗(𝐤⋅𝛔) dirac cone term
	Hᵣ = α*τ₁⊗(km[1]*σ₁ .+ km[2]*σ₂ .+ km[3]⊗σ₃)*√(abs(1-km[1]-km[2]))
	
	# Mass term
	#Hₘ = m*τ₃⊗I(2) 
	
	# Zeeman splitting + SoC cross term
	#Hᵦ = (β[1]⊗σ₁ .+ β[2]⊗σ₂ .+ β[3]⊗σ₃)
	Hᵦ = β[3]*I(2)⊗σ₃ .+ τ₃⊗(β[1]*σ₁.+β[2]*σ₂)
	
	H  = H₀ .+ Hᵦ .+ Hᵣ
	return H
end


klist = ["Γ","M","X₁","M","X₃","M","X₂","Γ"]
klist = ["Γ","M","X₁","Γ"]
C = π
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
plot3DSpectra(H,kdict["M"]+[0;0;0.00*π/a],a,160,160,0.1,0.1)
#plot3DSpectra(H,kdict["M"],a,160,160,0.1,0.1)
println("Done! Press ctrl+d to quit")

