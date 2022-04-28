

module ConstructHamiltonian
push!(LOAD_PATH, "./src/")
using LinearAlgebra
using Operators
using nearestNeighbors




function Hgen(Hstatic,edgeNNs,λ)
	R = Rvals(λ)
	NNs = nn1(λ,Rvals,t) #returns 1st nn hopping structs
	#one could modify NNs here
	H_hop, edgeNNs = nnHoppingMat(NNs,λ)
	Hstatic = H_hop # add on Hpot, Hcharge, H+U, etc here
	function H(k)
		N=2*λ^2
		H_edge = zeros(N,N)
		for NN in edgeNNs
			r = NN.r
			if(r[1]+r[2]>0)
				s = 1
			else
				s = -1
			end
			H_edge[NN.a,NN.b] += NN.t*exp(s*im*k⋅NN.r)
			#H_edge[NN.b,NN.a] += NN.t*exp(-im*k⋅NN.r)
		end
		return (H_edge.+Hstatic)
	end
	return H
end

end
#function hoppingModification(
