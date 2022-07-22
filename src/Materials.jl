#=module Materials
using LinearAlgebra
using Constants
using Operators
using UsefulFunctions
using SparseArrays

export Hopping, hoppingDict=#

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
	desc::String
end

function metalHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
	iorb = ia[5]
	t = (3*p.t)*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; ib = Int.(ia + di)
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        t = -1/2*p.t*I(2)
			pushHopping!(NNs, t, ia, ib, p)
                        #t = (p.ϵ₁ + 2*p.t)*(I(2))
                        #pushHopping!(NNs, t, ia, ia, p)
		end
	end
end

function insHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
	iorb = ia[5]
	t = (p.ϵ₁ + 3*p.t)*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; ib = Int.(ia + di)
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        t = -1/2*p.t*I(2)
			pushHopping!(NNs, t, ia, ib, p)
                        #t = (p.ϵ₁ + 2*p.t)*(I(2))
                        #pushHopping!(NNs, t, ia, ia, p)
		end
	end
end

function nextsite(iorb::Int)
    return (1-2*iorb)
end


function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0],""))
	push!(NNs,NN)
end


function weylHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        ib = ia; 
        #ib[5] += nextsite(iorb)
        #t = 3*p.t*(I(2))
	#pushHopping!(NNs, t, ia, ib , p)
        # for H₂ = τ₃⊗σ₀
        t = 3*nextsite(iorb)*p.t*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; 
                        # for Hweyl = τ₁⊗k⋅σ
                        di[5] = nextsite(iorb); 
                        ib = Int.(ia + di)
                        t = (-im/2)*dir*p.t*σ[ax]
			#Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			#δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
                        #t = nextsite(iorb)*(-im/2)*dir*p.t*σ[ax]
			pushHopping!(NNs, t, ia, ib, p)
			# for normal hopping term in hamiltonian
			ib[5] = iorb; 
			
			t = -(1/2)*nextsite(iorb)*p.t*(I(2))
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
end


hoppingDict = Dict{String,Function}(
                                    "mtjweyl"=>weylHopping,
                                    "weyl"=>weylHopping,
                                    "wins"=>insHopping,
                                    "ins"=>insHopping,
                                    "insulator"=>insHopping,
                                    "metal"=>metalHopping
                                   )


#end
