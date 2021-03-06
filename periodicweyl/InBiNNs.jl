push!(LOAD_PATH, "../src/")

module InBiNNs

using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
#using MaterialParameters
#using Materials
using Bands

export nnHoppingMat, genNNs, pruneHoppings, RvalsGen

# Takes in the parameters p, index vector (ix,iy,iz,isite (in unit cell), and iorb)
# returns the site-index in the full hamiltonian
function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	#ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4]; iorb = ivec[5]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	#ix = mod(ivec[1],p.nx); iy = mod(ivec[2],p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

include("../src/Materials.jl")
# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
function xyztor(p,ivec)
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4];
	δdict = Dict(0 => p.A*[0.0; 0.0; 0.0], #In 1 
		     1 => p.A*[0.5; 0.5; 0.5]) #In 2
		     #2 => p.A*[0.0; 0.5; 0.8975-0.5], #Bi 1
		     #3 => p.A*[0.5; 0.0; 1.10248-0.5]) #Bi 2
	δ = δdict[isite]
	R = p.a₁*ix + p.a₂*iy + p.a₃*iz + δ
        #println("ivec = $ivec, Rpos = $R")
        return R
end

function RvalsGen(p)
	N = p.n*p.nsite
	R = Vector{Vector{Float64}}(undef,N)
	for ix = 0:(p.nx-1)
		for iy = 0:(p.ny-1)
			for iz = 0:(p.nz-1)
				for isite = 0:(p.nsite-1)
					iR = Int(1 + isite + ix*p.nsite + iy*p.nx*p.nsite + iz*p.ny*p.nx*p.nsite)
					#println("$ix $iy $iz $isite iR $iR")
					ivec = Int.([ix,iy,iz,isite])
                                        Rval = xyztor(p,ivec)
                                        R[iR] = deepcopy(Rval)
                                        #if(Rval[3] > p.a₃[3]*(p.nz-1))
                                        #    println("ivec = $ivec, Rpos = $Rval")
                                        #end
                                end
			end
		end
	end
        #println("maximum R = $(maximum([pos[3] for pos in R])), expected height = $((p.nz-1)*p.a₃[3])")
	return R # needs ⊗I(p.norb)⊗I(2) for full (spinful) hilbert space
end

# Defines a cᵦ†cₐ term 
#=mutable struct Hopping
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
end=#


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
			#println("NN hopping = $(NN.t)")
			#println("indices a,b = $(NN.a), $(NN.b)")
			#show(NN)
			H[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	#println("H₀ = $H")
	return H, edgeNNs
end

# Add a bond to the list of bonds, given some list of bonds, coefficient in spin basis, index of both sites, and param list
#=function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0],""))
	#show(NN)
	if(abs(ib⋅[0,0,1,0,0]) ≈ 1)
		dir = "z"
	elseif(abs(ib⋅[0,1,0,0,0]) ≈ 1)
		dir = "y"
	elseif(abs(ib⋅[1,0,0,0,0]) ≈ 1)
		dir = "x"
	else
		dir = ""
	end
	if((ib⋅[1,1,1,0,0]) < 0)
		sgn = " in -"
	elseif(ib⋅[1,1,1,0,0] ≈ 0)
		sgn = ""
	else
		sgn = " in +"
	end
	if(t[1,1] ≈ t[2,2])
		if(t[1,1] ≈ 0)
			if(t[1,2] ≈ -t[2,1])
				σ = "σ₂"
				C = t[2,1]/im
			elseif(t[1,2] ≈ t[2,1])
				σ = "σ₁"
				C = t[1,2]/1
			else
				σ = "unknown"
				C = "unknown"
			end
		elseif(t[1,2] ≈ 0)
			σ = "σ₀"
			C = t[1,1]
		else
			σ = "unknown"
			C = "unknown"
		end
	elseif(t[1,1] ≈ -t[2,2])
		σ = "σ₃"
		C = t[1,1]
	else
		σ = "unknown"
		C = "unknown"
	end
	desc = "ia = $(NN.ia) to ib = $(NN.ib)"
	#desc = "Hop as τ = |$(ib[5])$sgn$dir><$(ia[5])| with t = $(round(C,sigdigits=3))$σ, |$(NN.b)><$(NN.a)|"
	NN.desc = desc
	#println(desc)
	push!(NNs,NN)
end=#

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


function genNNs(p) # all of the terms in the hamiltonian get added here, get back the relevant bonds
	n = p.n
	NNs = Hopping[]
	
	function nextsite(isite::Int)
		return -2*isite + 1
		#if == 0 return 1
		#if == 1 return -1
		#return mod(isite+1,p.nsite)
	end
	orb = 0;
	# loop over each unit cell site in the superlattice
	for iy = 0:(p.ny-1)
		for iz = 0:(p.nz-1)
			for ix = 0:(p.nx-1)
				isite = 0
				for iorb = 0:(p.norb-1)
					ia = copy([ix,iy,iz,isite,iorb]);
					#pushHopping!(NNs, -nextsite(iorb)*3*p.t*I(2), ia, ia, p)
                                        hoppingType! = hoppingDict[p.deviceMaterial]
                                        hoppingType!(p,NNs,ia)
                                        #=
                                        t = 3*nextsite(iorb)*p.t*(I(2))
					pushHopping!(NNs, t, ia, ia, p)
					for ax = 1:3
						for dir = [-1,1]
							# for weyl term in hamiltonian
							di = zeros(5); di[ax] = dir; di[5] = nextsite(iorb); ib = Int.(ia + di)
							Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
							δ = Rb - Ra
							# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
							#println("δR$ax = $dR")
							#pushHopping!(NNs, (p.vf*-im*ħ/q)*(1/(2*dR))⊗I(2), ia, ib, p)
							#t = zeros(ComplexF64,2,2)
							#for ax = 1:3
								#dR = abs(Rb[ax] - Ra[ax])
								#if(dR ≈ 0)
								#	println("Uh oh, some distances are 0...")
								#	println("di = $di")
								#end
							#t = (1/(2*δ[ax]))*σ[ax]
							t = (-im/2)*dir*p.t*σ[ax]
							#t = -im*p.t*(δ[ax]/(2*norm(δ)))*σ[ax]
							#t = -im*(p.vf*ħ/q)*(1/(2*δ[ax]))*σ[ax]
							if(any(isnan,t)==true)
								throw(DomainError(t, "Something broken in hamiltonian definition! Returning NaN"))
								return
							end
							pushHopping!(NNs, t, ia, ib, p)
							# for normal hopping term in hamiltonian
							ib[5] = iorb; 
							
							t = -(1/2)*nextsite(iorb)*p.t*(I(2))
							pushHopping!(NNs, t, ia, ib, p)
						end
					end=#
				end
			end
		end
	end
	# now fix the designation for the vectors that hop out of the lattice
	# Will connect them around later using bloch's theorem to generate H(k) function
	for NN in NNs
		#println("pre hop ($(NN.ia) to $(NN.ib)) = $(NN.a) to $(NN.b)")
		ib = [NN.ib[1],NN.ib[2],NN.ib[3]]
		# Δ(ib,ib reflected back into 1st lattice)
		pib = ib - [mod(ib[1],p.nx),mod(ib[2],p.ny),mod(ib[3],p.nz)]
		if(pib⋅pib != 0) # if vector is distinctly outside of 1st lattice
			NN.N = Int.([round(pib[1]/(p.nx)),round(pib[2]/p.ny),round(pib[3]/p.nz)])
			NN.b = xyztoi(p,NN.ib, NN.N)
			NN.edge = true
			#println("$(NN.N)")
		end
		#println("NN hop $(NN.a) to $(NN.b) = $(NN.desc)")
		#println("edge? = $(NN.edge)")
	end
	return NNs
end


function pruneHoppings(NNs, type)
	if("x" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[1]!=0,NNs))
	end
	if("y" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[2]!=0,NNs))
	end
	if("z" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[3]!=0,NNs))
	end
	return NNs
end

end
