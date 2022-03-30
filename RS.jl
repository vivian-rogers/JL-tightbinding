push!(LOAD_PATH, "./src/")


using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using MaterialParameters
using Bands
using PyPlot



# band calculation parameters
klist = ["W","L","Γ","X","W","K"] 	#passes them to Bands, which matches them with numerical values
				#double-check values of kpoints in band plotting routine in src/Bands.jl if bands look weird
nk = 200 			#number of kpoints between given high-sym points
proj = true #project eigenstate onto given basis vector? might be broken, disregard
projKet = [1;0]⊗[0; 0; 1] #projects onto orbital |a>⊗|↑> + |a>⊗|↓>

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
	
	nₐ = 2; nₒ = 3; nₛ = 1 #number of atoms, orbitals, and spins
	# basis is [|Sc> |N>]⊗[|s> |p> |d>]⊗[|↑> |↓>]
	n = nₐ*nₒ*nₛ # number of bands considered
	H = zeros(ComplexF64, n,n)
 	
	# Hamiltonian  	of a single Sc atom, in basis |orbital a> |orbital b>
	# Modify as needed
	#e1 = 1; e2 = 1.5; e3 = 1.2; e4 = 3
	εSc = eV*Diagonal([0, 0., 0.9]); εN = eV*Diagonal([-2,0,10]);
	ε = [1 0; 0 0]⊗εSc + [0 0; 0 1]⊗εN; 
	#ε = eV*Diagonal([0.9 0]); #in basis of Sc, N
	t_nn_sp = 0.5*eV;
	t_nn_sd = 0.7*eV;
	t_nnn_dd = 0.1*eV;
	t_nnn_ss = 0.005*eV;
	t_sc_pd = 0.01*eV;
	t_sc_sd = 0.4*eV
	t_sc_sp = 0.8*eV

	τSc = eV*[
	       	0	t_sc_sp	t_sc_sd
	       	t_sc_sp	0	t_sc_pd
	       	t_sc_sd	t_sc_pd	0
	       ]
	τN = eV*[
		0	0	0
		0	0	0
		0	0	0
		]
	τScN = eV*[
			0	t_nn_sp	t_nn_sd
			t_nn_sp	0	0
			t_nn_sd	0	0
		]
	#t2 = 0.3*eV;
	
	#orbital hopping matrix between NNs
	τ = Diagonal([1,0])⊗τSc .+ Diagonal([0,1])⊗τN;
	
	#atomic energies + intra-unitcell hopping
	H .= ε .+ τ; #for within the unit cell
	
	d1Vects = nn1_vects(a) #gives list of 1st nearest neighbor peaks
	d2Vects = nn2_vects(a) #gives list of 2nd nearest neighbor peaks
	for r in d1Vects
		H .+= exp(im*k⋅r)*(τ₁⊗τScN)
	end
	for r in d2Vects
		H .+= exp(im*k⋅r)*(Diagonal([1,0])⊗(Diagonal([t_nnn_ss, 0, t_nnn_dd])))
	end
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
#show((1)*i for i in nn1_vects(a))
#show(nn2_vects(a))
#show(nn2_vects)
E, Estates = getBands(klist, kdict, nk, a, H)
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

