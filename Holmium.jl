push!(LOAD_PATH, "./src/")


using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using HolmiumParameters
using Bands
using PyPlot



# band calculation parameters
klist = ["Γ","M","K","Γ"] 	#passes them to Bands, which matches them with numerical values
				#double-check values of kpoints in band plotting routine in src/Bands.jl if bands look weird
nk = 400 			#number of kpoints between given high-sym points
proj = false #project eigenstate onto given basis vector? might be broken, disregard
projKet = [1;0]⊗[1; 0]⊗[1;0] #projects onto orbital |a>⊗|↑> + |a>⊗|↓>

function H(k)
	
	# Hamiltonian(k) generator
	# LaTeX subscripts can be typed in vim, Ex: \_s + tab -> ₛ
	# without julia/vim editor, simply copy/paste characters or delete them and replace with ASCII
	# typing \otimes + tab will give kronecker product, ⊗
	# typing I(n) will give n x n identity matrix 
	# typing \sigma + tab and \Sigma + tab will give σ and Σ
	#
	# Values such as hopping parameters t, γ, β, global constants, etc, can be stored in src/Constants.jl and src/HolmiumParameters.jl
	# functions such as the rotate-around-z function exist in src/UsefulFunctions.jl
	# .+, .=, .+=, etc will access arrays and matrices in traditional way, use this generally
	
	nₐ = 2; nₒ = 2; nₛ = 2 #number of atoms, orbitals, and spins
	n = nₐ*nₒ*nₛ # number of bands considered
	H = zeros(ComplexF64, n,n)
 	
	# Hamiltonian of a single Ho atom, in basis |orbital a> |orbital b>
	# Modify as needed	
	εₐ = [
		1 0; 
		0 1.5
	     ]*eV 

	#orbital hopping matrix between NNs
	tₕ = [
	      t*0.3 	t*0.1; 
	      t*0.1 	t*1.2
	      ]
	
	#atomic energies + intra-unitcell hopping
	#H .= I(nₐ)⊗εₐ⊗I(nₛ)
	H .= (I(nₐ)⊗εₐ .+ [0 1; 1 0]⊗tₕ)⊗I(nₛ)

	# 2D hopping 
	for i = 0:2
		# r₁ vector is defined in src/HolmiumParameters, may warrant reformulation in terms of lattice vectors	
		r = rotate(i*π/3)*r₁
		H .+= exp(im*k⋅r)*I(nₐ)⊗tₕ⊗I(nₛ)
	end

	#out-of-plane hopping
	#H .+= exp(im*k⋅r₂)*[0 1; 0 0]⊗tₕ⊗I(nₛ)
	#H .+= exp(-im*k⋅r₂)*[0 0; 1 0]⊗tₕ⊗I(nₛ)
	#=for i = 0:2
		# |1>,|2>; 
		r = rotate(i*2*π/3)*r₂
		H .+= exp(im*k⋅r)*[0 1; 0 0]⊗tₕ⊗I(nₛ)
		r = rotate(i*2*π/3)*r₃
		H .+= exp(im*k⋅r)*[0 0; 1 0]⊗tₕ⊗I(nₛ)
	end =#
	
	#example zeeman term in Z
	H .+= 1*eV*I(nₐ)⊗I(nₒ)⊗σ₃
	return H
end


#println("Volume of unit cell:")
#A = [
#	a cosd(60)*a 0;
#	0 sind(60)*a 0;
#	0 0          c
#    ]
#display(A)

#println("\nVolume = $(det(A)) m^3")
println("Getting eigenvalues of holmium between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, nk, a, H)
if(proj)
	projStates = project(projKet, Estates)
end


if(proj)
	println("Plotting with projection onto given basis ket...")
	plotBands(klist,nk,E, projStates)
else
	println("Plotting...")
	plotBands(klist,nk,E)
end
println("Done! Press ctrl+d to quit")

