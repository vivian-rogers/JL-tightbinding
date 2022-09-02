module PWHamiltonian

using LinearAlgebra
using SparseArrays
using UsefulFunctions
using Constants
using Operators
using Distributed

export gGrid, gridOffset, Hβgen, ConstructHamiltonian


mutable struct Hopping
	a::Int # orbital/site index 1
	b::Int # orbital/site index 2 with PBC
	ia # index vector of site A
	ib # index vector of site B without PBC
	t  # hopping parameter affiliated with c†₂c₁ in spin basis. (i.e., t*I(2) or t*σ₊ might make sense)
	desc::String
end

function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = gtoi(p,ia); b = gtoi(p,ib);
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,t,""))
	push!(NNs,NN)
end

function gtoi(p,gvec)
        # gives index in matrix for g vector (-gmax:gmax)
        iG = gvec[1]; iorb = gvec[2]
        return iG + iorb*p.nG
end

function ConstructHamiltonian(p,M)
    Hᵦ = Hβgen(p,M)
    Hfree = freeElectron(p)
    Hweyl = weylH(p)
    H₀ = (10^-4 + 10^-6*im)*I(p.nG*p.norb*2) .+ Diagonal(p.μ_disorder*(rand(p.nG*p.norb*2).-0.5))
    function H(k::Union{Vector{ComplexF64},Vector{Float64}})
        k = k.+[1;1;1]*10^-8 # offset k slightly to avoid weird stuff
        H =  dropzeros(Hweyl(k) .+ Hᵦ .+ H₀)
        #H =  Hweyl(k) .+ I(p.ng*p.norb)⊗σ₂
        #=println("Hfree")
        display(Hfree(k))
        println("H periodic magnetic")
        display(Hᵦ)=#
        #H =  Hfree(k) .+ I(p.nG*p.norb)⊗σ[3]
        #H =  Hfree(k) .+ Hᵦ .+ H₀
        #H = Hβ(k) .+ Hfree(k) .+ 2*Diagonal(ones(p.ng*2*p.norb))
        return H
    end
end

mutable struct Gvec
    ng::Vector{Int}
    G::Vector{Float64}
    E::Float64
    i::Int # index in hamiltonian
end

function gGrid(p::NamedTuple)
    Gvecs, nG = gGrid(p.B,p.maxG,p.gcut,p.vf)
    return Gvecs
end

function gGrid(B::Matrix, maxG::Vector{Int} = [10,10,10], gcut = 50.0*eV, vf = 10^6)
        gvecs = Gvec[]
        #println("Generating $(2*maxG + [1,1,1]) G grid...")
        for gx = -maxG[1]:maxG[1]
            for gy = -maxG[2]:maxG[2]
                for gz = -maxG[3]:maxG[3]
                    ng = [gx;gy;gz]
                    G = B*ng
                    #E = vf*(ħ/q)*√(G⋅G)
                    E = ħ^2*(G⋅G)/(2*m₀*q)
                    #E = ħ^2*(G⋅G)/(2*m₀*q)
                    push!(gvecs,deepcopy(Gvec(ng,G,E,0)))
                end
            end
        end
        # remove all g vectors above a threshold
        if(typeof(gcut)==Float64)
            deleteat!(gvecs, findall(g->(g.E > gcut),gvecs))
        else
            deleteat!(gvecs, findall(g->(norm(g.G) > norm(B*gcut)),gvecs))
        end
        #println("Using $(size(gvecs)[1]) / $(prod(maxG*2 + [1,1,1])) G vectors")
        # now add their indices
        for ig in eachindex(gvecs)
            gvecs[ig].i = deepcopy(ig)
        end
        nG = maximum([g.i for g in gvecs])
        return gvecs, nG
end




function gridOffset(p::NamedTuple, ivec::Vector{Int})
    i = (p.maxG+[1,1,1]) + ivec
    ind = CartesianIndex(i[1],i[2],i[3])
    #println("Gvec # $(ivec), index # $(ind)")
    return ind
end

#=function weylH2(p)
    gvecs = gGrid(p)
    function Hweyl(k::Vector{Float64})
   
end=#


function weylH(p)
    gvecs = gGrid(p)
    function Hweyl(k::Union{Vector{ComplexF64},Vector{Float64}})
        Hw = spzeros(ComplexF64,p.ng*p.norb*2,p.ng*p.norb*2)
        # build up a tridiagonal matrix corresponding to the diagonals in the weyl hamiltonian
        # in hilbert space |g>⊗|spin>
        
        # linear dirac hamiltonian
        #=NNs = Hopping[]
        for iG in eachindex(gvecs)
            G = gvecs[iG]
            for iorb = 0:(p.norb-1)
                ia = [G.i, iorb]
                ib = [G.i, Int(!(Bool(iorb)))] # maps 0->1, 1->0 for τₓ orbital
                a = gtoi(p,ia); b = gtoi(p,ib)
                t = p.t*matdot(sin.(k+G.G),(σ))
                H[2*b-1, 2*a-1] = t[1,1]
                H[2*b  , 2*a-1] = t[2,1]
                H[2*b-1, 2*a  ] = t[1,2]
                H[2*b  , 2*a  ] = t[2,2]
            end
        end=#
        NNs = Hopping[]
        for iG in eachindex(gvecs)
            G = gvecs[iG]
            for iorb = 0:(p.norb-1)
                ia = [G.i, iorb]
                ib = [G.i, Int(!(Bool(iorb)))] # maps 0->1, 1->0 for τₓ orbital
                a = gtoi(p,ia); b = gtoi(p,ib)
                t = p.t*matdot(sin.(p.A*(k+G.G)),(σ))
                Hw[2*b-1, 2*a-1] = t[1,1]
                Hw[2*b  , 2*a-1] = t[2,1]
                Hw[2*b-1, 2*a  ] = t[1,2]
                Hw[2*b  , 2*a  ] = t[2,2]
            end
        end
        # now add the mass term 
        Hm = spzeros(ComplexF64,p.ng*p.norb,p.ng*p.norb)
        NNs = Hopping[]
        for iG in eachindex(gvecs)
            G = gvecs[iG]
            for iorb = 0:(p.norb-1)
                ia = [G.i, iorb]
                ib = [G.i, iorb] # maps 0->1, 1->0 for τₓ orbital
                a = gtoi(p,ia); b = gtoi(p,ib)
                t = ((-1)^iorb)*p.t*(3 - cos.(p.A*(k+G.G))⋅[1;1;1])
                Hm[b, a] = t
            end
        end
        return (Hw .+ p.m*Hm⊗I(2))
    end
    return Hweyl
end


function freeElectron(p)
    gvecs = gGrid(p)
    function Hfree(k::Vector{Float64})
        H = zeros(p.ng)
        for g in gvecs
            E = ħ^2*norm(g.G + k)^2/(2*m₀*q)
            H[g.i] = deepcopy(E)
        end
        return Diagonal(H)⊗I(p.norb)⊗I(2)
    end
    return Hfree
end

function HConstruct(NNs::Vector{Hopping},p::NamedTuple)
	N = p.nG*p.norb
	H = spzeros(ComplexF64,2*N,2*N)
	edgeNNs = Hopping[]
	for NN in NNs
            # at site 1->2 would be a_2,1
            # with site ⊗ spin, would be 
            #H[NN.a,NN.b] = NN.t
            #println("NN hopping = $(NN.t)")
            #println("indices a,b = $(NN.a), $(NN.b)")
            #show(NN)
            H[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
            H[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
            H[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
            H[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
	end
	#println("H₀ = $H")
	return H
end


function Hβgen(p::NamedTuple,M)
    gvecs = gGrid(p)
    Hb = spzeros(ComplexF64, p.ng*p.norb*2, p.ng*p.norb*2)
    Mx = M[1]; My = M[2]; Mz = M[3]
    NNs = Hopping[]
    for G in gvecs
        for G2 in gvecs #now do the convolution
            for iorb = 0:(p.norb-1)
                ia = [G2.i,iorb]
                t = zeros(ComplexF64,2,2)
                t .+= M[1][gridOffset(p,G.ng-G2.ng)]*I(2)
                for ax = 2:4
                    Mᵢ = M[ax]
                    t .+= Mᵢ[gridOffset(p,G.ng-G2.ng)]*σ[ax-1]
                end
                ib = [G.i,iorb]
                pushHopping!(NNs, t, ia, ib, p) 
            end
        end
    end
    Hb = dropzeros(HConstruct(NNs,p))
    return Hb
end


end
