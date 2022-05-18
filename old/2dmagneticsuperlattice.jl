push!(LOAD_PATH, "../src/")

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
	H₀ = zeros(ComplexF64,4,4)
	# add interactions of form |A><B| hopping out of unit cell
	# change sign of hopping parameter to be + (-) in +i (-i) hopping, add tensor product with right hopping direction
	for di = [
		 [1,0,1],[0,1,1],[0,0,1],[1,1,1],
		 #[1,0,0],[0,1,0],	[1,1,0]
		[1,0,0],[0,1,0],[0,0,0],[1,1,0]
		 ]
		for ax = 1:3
			t = di[ax]*2 - 1
			H₀ .+= exp(k⋅di)*[0 t; 0 0]⊗σ[ax]
		end
	end
	# add interactions of form |B><A| hopping out of unit cell
	for di = [
		 [0,0,0],[0,-1,0],[-1,-1,0],[-1,0,0],
		# 	[0,-1,0],[-1,-1,0],[-1,0,0],
		[0,0,-1],[0,-1,-1],[-1,-1,-1],[-1,-1,-1]
		 ]
		for ax = 1:3
			t = di[ax]*2 + 1
			H₀ .+= exp(k⋅di)*[0 0; t 0]⊗σ[ax]
		end
	end
	H = H₀	
	return Hermitian(H)
end

#=function H(k) fcc weyl
	#in basis |1> |2> for the 2 atom unit cell
	m = 0*eV # mass term at M point
	t = 1.0*eV # hopping parameter
	#t₂ = 0.4*eV
	#t₃ = 0.2*eV
	#α = 0.01*eV
	#α₂ = 0.0*eV
	ε = 0
	r = [1 0 0]
	kx = k[1]; ky = k[2]; kz = k[3];
	H₀ = zeros(ComplexF64,4,4)
	#H₀ .+= [0 1; -1 0]⊗(σ₁.+σ₂.+σ₃)
	#H₀ .+= [0I(2) σ₁.+σ₂.+σ₃; -(σ₁.+σ₂.+σ₃) 0]
	# connections from site 1
	for di = [
		 [1,0,1],[0,1,1],[0,0,1],[1,1,1],
		 #[1,0,0],[0,1,0],	[1,1,0]
		[1,0,0],[0,1,0],[0,0,0],[1,1,0]
		 ]
		for ax = 1:3
			t = di[ax]*2 - 1
			H₀ .+= exp(k⋅di)*[0 t; 0 0]⊗σ[ax]
		end
	end
	for di = [
		 [0,0,0],[0,-1,0],[-1,-1,0],[-1,0,0],
		# 	[0,-1,0],[-1,-1,0],[-1,0,0],
		[0,0,-1],[0,-1,-1],[-1,-1,-1],[-1,-1,-1]
		 ]
		for ax = 1:3
			t = di[ax]*2 + 1
			H₀ .+= exp(k⋅di)*[0 0; t 0]⊗σ[ax]
		end
	end
	H = H₀	
	return Hermitian(H)
end=#


klist = ["Γ","M","X₁","M","X₃","M","X₂","Γ"]
klist = ["Γ","M","X₁","Γ"]
C = 2π
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

