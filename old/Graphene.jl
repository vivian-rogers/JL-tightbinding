push!(LOAD_PATH, "./src/")


using Plots
using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using GrapheneParameters
using Bands



⊗(A,B) = kron(A,B)


function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	
	H₂₁ = [0 0
	       t 0] #hop from A to B
	
	H₁₂ = [0 t
	       0 0] #hop from B to A
	
	H₀ = ε*I(2) .+ H₂₁ .+ H₁₂

	# maybe add in peierls phase sub for 𝐀 = [0,Bx,0]?
	H = H₀ .+ exp(im*k⋅r₁)*H₂₁ .+ exp(im*k⋅r₂)*H₂₁ .+ exp(-im*k⋅r₁)*H₁₂ .+ exp(-im*k⋅r₂)*H₁₂ 
	
	β = 0.0*eV
	H = H⊗I(2) .+ β*I(2)⊗σ₃
	return H
end


klist = ["Γ", "M", "K", "Γ"]

C = 4*π/a
kdict = Dict(
	    "Γ" => C*[0;   0;   0],
	    "K" => C*[1/2; tand(30)/2; 0],
	    "M" => C*[1/2; 0;   0],
	    "M-" => C*[-1/2; 0;   0],
	    "K'" =>C*[1/2; -tand(30); 0],	
	    "H" => C*[1/3;1/3;1/2]
	    )	    
nk = 1028
println("Getting eigenvalues of graphene between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, kdict, nk, a, H)
#display(27.2*E)
println("Plotting...")
plotBands(klist,nk,E)
println("Done! Press ctrl+d to quit")

