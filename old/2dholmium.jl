push!(LOAD_PATH, "./src/")


using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using HolmiumParameters
using Bands

function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	H₀ = [ε] #hamiltonian of the unit cell
	H₁ = [t] #so <2|1> has a hopping energy
	#H = H₀ .+ exp(im*k⋅(R(0)*a₁))*H₁ .+ exp(im*k⋅(R(a₂))*H₁ .+ exp(im*k⋅(a₁+a₂))*H₁ .+
	#	exp(im*k⋅(R(180)*a₁))*H₁ .+ exp(-im*k⋅a₂)*H₁ .+ exp(-im*k⋅(a₁+a₂))*H₁ 	
	H = H₀
	for i in 0:2
		r = 2/3*rotate(i*π/3)*a₁
		H += exp(im*k⋅r)*H₁
	end
	γ = 0.3*eV
	Hₛₒ = γ*(sin(π*k[1])*σ₃ .+ sin(π*k[2])*σ₁)
	display(Hₛₒ .+ Hₛₒ')
	H = H⊗I(2) .+ I(1)⊗(Hₛₒ .+ Hₛₒ')
	return H
end


klist = ["Γ", "M", "K", "Γ"]
nk = 1028
println("Getting eigenvalues of graphene between k = ")
show(klist)
println("...")
E = getBands(klist, nk, a, H)
#display(27.2*E)
println("Plotting...")
plotBands(klist,nk,E)
println("Done! Press ctrl+d to quit")

