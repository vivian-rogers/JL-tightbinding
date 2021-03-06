module GenNNs

using MaterialParameters

export NNatoms

struct AtomNN
	A1 # atom 1
	A2 # atom 2 
	δ # radius from A1 to A2
	N # unit cell of second atom, [Nx, Ny, Nz]⋅[a₁, a₂, a₃]
end

function NNatoms()
	nn1_dist = a/2;
	#AtomNNs = AtomNN[]
	AtomNNs = Array{AtomNN, 1}(undef, (6+8)*4) # each atom will have 6 nn₂s, 8 nn₁s
	n = 1
	for atom = 1:4
		for ax = 1:3
			for dir = [-1,1]
				A1 = Atoms[atom]
				A2 = Atoms[mod((atom+2-1),4)+1]
				δ = zeros(3)
				δ[ax] = dir*(a/2)
				#N = zeros(3)
				
				# lot going on here. This gives distance from imagined periodic
				# atom to the base atom in the unit cell. Take this distance and 
				# multiply by [a1,a2,a3]^-1 (so B) but include (2*π) because
				# this isn't the reciprocal lattice
			
				N = Int.(round.((2*π)^(-1)*B*(δ + A1.R - A2.R)))
				AtomNNs[n] = AtomNN(A1, A2, δ, N)
				n += 1
			end
		end
		for d1 = [-1,1]
			for d2 = [-1,1]
				for d3 = [-1,1]
					dn = d1*d2*d3
					δ = (a/4)*[d1;d2;d3]
					A1 = Atoms[atom]
					A2 = Atoms[mod(atom+dn-1,4)+1]
					N = Int.(round.((2*π)^(-1)*B*(δ + A1.R - A2.R)))
					AtomNNs[n] = AtomNN(A1, A2, δ, N)
					n += 1	
				end
			end
		end
	end
	return AtomNNs
end

end
