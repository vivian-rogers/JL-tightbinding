push!(LOAD_PATH, "../src/")

module InBiNNs

using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
#using MaterialParameters
using Bands


export nnHoppingMat, genNNs

# Takes in the parameters p, index vector (ix,iy,iz,isite (in unit cell), and iorb)
# returns the site-index in the full hamiltonian
function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	#ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4]; iorb = ivec[5]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
function xyztor(p,ivec)
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4];
	δdict = Dict(0 => p.A*[0.0; 0.0; 0.0]) #In 1 
		     #1 => p.A*[0.5; 0.5; 0.0]) #In 2
		     #2 => p.A*[0.0; 0.5; 0.8975-0.5], #Bi 1
		     #3 => p.A*[0.5; 0.0; 1.10248-0.5]) #Bi 2
	δ = δdict[isite]
	return p.a₁*ix + p.a₂*iy + p.a₃*iz + δ
end

# Defines a cᵦ†cₐ term 
mutable struct Hopping
	a::Int # orbital/site index 1
	b::Int # orbital/site index 2 with PBC
	ia # index vector of site A
	ib # index vector of site B without PBC
	ra # location of atom a
	rb # location of atom b
	r # radius from a to b
	t  # hopping parameter affiliated with c†₂c₁ in spin basis. (i.e., t*I(2) or t*σ₊ might make sense)
	edge::Bool # does this hop off the edge of the superlattice?
	N # vector describing the [n₁;n₂;n₃]⋅[a₁;a₂;a₃] superlattice unit cell of site ib
end


# Generate H₀ and make a list of edge bonds for generating H(k)
function nnHoppingMat(NNs,p)
	N = p.n*p.nsite*p.norb
	H = zeros(ComplexF64,2*N,2*N)
	edgeNNs = Any[]
	for NN in NNs
		if(NN.edge == true)
			push!(edgeNNs,deepcopy(NN))
		else
			# at site 1->2 would be a_2,1
			# with site ⊗ spin, would be 
			#H[NN.a,NN.b] = NN.t
			H[2*NN.b-1, 2*NN.a-1] = NN.t[1,1]
			H[2*NN.b  , 2*NN.a-1] = NN.t[2,1]
			H[2*NN.b-1, 2*NN.a  ] = NN.t[1,2]
			H[2*NN.b  , 2*NN.a  ] = NN.t[2,2]
		end
	end
	return H, edgeNNs
end

# Add a bond to the list of bonds, given some list of bonds, coefficient in spin basis, index of both sites, and param list
function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0])
	push!(NNs,NN)
end

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


function genNNs(p) # all of the terms in the hamiltonian get added here, get back the relevant bonds
	n = p.n
	NNs = Hopping[]
	
	In1 = 0; In2 = 1; Bi1 = 2; Bi2 = 3; px = 0; py = 1
	# loop over each unit cell site in the superlattice
	for iy = 0:(p.ny-1)
		for iz = 0:(p.nz-1)
			for ix = 0:(p.nx-1)
				for isite = 0:(p.nsite-1)
					# loop over each site in the given ix,iy,iz unit cell and generate cᵦ^†cₐ terms on hamiltonian
					# note: Δi represents the offset to the index vector, the hopping term loop is defined using such
					# key: di = [Δi_x, Δi_y, Δi_z, iᵦ - iₐ (dif of atom indices in unit cell), Δi_orbital
					# Then ia and ib represent hamiltonian matrix index for hilbert space:
					# [Unit cell] ⊗ [site] ⊗ [orbital px or py]
					# Hamiltonian is now spinless, see p.t₂*I(2) in line 108. Can use p.t₁*σ₁, say. 
					ia = [ix,iy,iz,isite,0];
					for ax = [1,2,3]
						for dir = [-1,1]
							ib = ia; di = zeros(5); di[ax] = dir; ib += di; ib = Int.(ib)
							Ra = xyztor(p,ia); Rb = xyztor(p,ib); dR = dir*norm(Rb - Ra) 
							# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
							pushHopping!(NNs, (p.vf*-im*ħ/q)*dir*(1/(2*dR))⊗σ[ax], ia, ib, p)
						end
					end
				end
			end
		end
	end
	# now fix the designation for the vectors that hop out of the lattice
	# Will connect them around later using bloch's theorem to generate H(k) function
	for NN in NNs
		ib = [NN.ib[1],NN.ib[2],NN.ib[3]]
		# Δ(ib,ib reflected back into 1st lattice)
		pib = ib - [mod(ib[1],p.nx),mod(ib[2],p.ny),mod(ib[3],p.nz)]
		if(pib⋅pib != 0) # if vector is distinctly outside of 1st lattice
			NN.N = Int.([round(pib[1]/(p.nx)),round(pib[2]/p.ny),round(pib[3]/p.nz)])
			show(NN.N)
			NN.b = xyztoi(p,NN.ib, NN.N)
			NN.edge = true
			#println("$(NN.N)")
		end
	end
	return NNs
end


end
