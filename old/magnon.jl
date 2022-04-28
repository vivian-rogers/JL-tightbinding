push!(LOAD_PATH, "./src/")


using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using HolmiumParameters
using Bands


klist = ["Γ","M","K","Γ","A"]
nk = 2000
proj = true
#projKet = [1; 0] #projects onto orbital |a>⊗|↑> + |a>⊗|↓>
#projKet = [1;1]⊗[0; 1]⊗[1;0] #projects onto orbital |a>⊗|↑> + |a>⊗|↓>
projKet = [1;0]⊗[1; 1] #projects onto orbital |1>⊗|↑> + |a>⊗|↓>

function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	nOrbitals = 2; nAtoms = 2; nSpins = 1
	n = nOrbitals*nAtoms*nSpins # number of bands considered
	H = zeros(ComplexF64, n,n)
 	
	#hamiltonian of the Ho atom in basis |atom, orbital>
 	Hₐ = I(nAtoms)⊗[
		   0.06 0; 
		   0.07 0
		  ] 

	#orbital hopping matrix between NNs
	tₕ = [
	      t 	t*0.1; 
	      t*0.1 	t*2
	      ]
	
	#atomic energies + intra-unitcell hopping
	H .= Hₐ .+ [0 1; 1 0]⊗tₕ

	# 2D hopping 
	for i in 0:2
		#add back to same atom with bloch phase difference
		r = rotate(i*π/3)*r₁
		H .+= exp(im*k⋅r)*I(2)⊗tₕ
	end

	#out-of-plane hopping
	for i in 0:2
				# |1>,|2>; 
		r = rotate(i*2*π/3)*r₂
		H .+= exp(im*k⋅r)*[0 1; 0 0]⊗tₕ
		r = rotate(i*2*π/3)*r₃
		H .+= exp(-im*k⋅r)*[0 0; 1 0]⊗tₕ
	end
	#H = H⊗I(2) + 0.1*eV*I(4)⊗S₃
	return H
end



println("Getting eigenvalues of graphene between k = ")
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

