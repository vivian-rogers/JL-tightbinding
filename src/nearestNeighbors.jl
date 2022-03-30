

#push!(LOAD_PATH, "./src/")
module nearestNeighbors

using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using MaterialParameters
using Bands


export nnHoppingMat, NNfromAtoms 

mutable struct Hopping
	a::Int # orbital/site index 1
	b::Int # orbital/site index 2
	r # radius from a to b
	t  # hopping parameter affiliated with c†₂c₁
	edge::Bool
	N # vector describing the [n₁;n₂]⋅[a₁;a₂] superlattice unit cell
end

function nnHoppingMat(NNs,N)
	H = zeros(ComplexF64,N,N)
	edgeNNs = Any[]
	for NN in NNs
		if(NN.edge == true)
			push!(edgeNNs,deepcopy(NN))
		else
			H[NN.a,NN.b] = NN.t
		end
	end
	return H, edgeNNs
end

function NNfromAtoms(NNAtoms)
	n = size(NNAtoms)[1]
	NNs = Array{Hopping, 1}(undef, n)
	i = 1
	for NNatom in NNAtoms
		i1 = NNatom.A1.i
		i2 = NNatom.A2.i
		#show([i1,i2])
		if(NNatom.N == [0;0;0])
			edge = false
		else
			edge = true
		end
		δ = NNatom.δ
		
		# this needs fixing when importing from wannier90
		# whole subroutine dedicated to importing hoppings; this is estimate
		#t = 3.0*exp(-(1/3)*norm(δ)/Å)
		t = tdict[[NNatom.A1.type,NNatom.A2.type]]
		NNs[i] = Hopping(i1,i2,δ,t,edge,NNatom.N)
		i += 1
	end
	return NNs
end

end
