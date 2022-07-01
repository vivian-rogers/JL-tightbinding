module Electrodes
using LinearAlgebra
using Constants
using Operators
using UsefulFunctions
using SparseArrays
using VectorPotential

export Electrode, genΣₖs

mutable struct Electrode
	xrange::Vector{Int}
	yrange::Vector{Int}
	zrange::Vector{Int}
	n::Int
	connectfrom::String # direction connecting fThreads.@threadsrom (i.e, "-x","+x",etc)
	type::String # material that the electrode is made of
	A::Function # magnitude of the exchange field
end


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
 



function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	#ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4]; iorb = ivec[5]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	#ix = mod(ivec[1],p.nx); iy = mod(ivec[2],p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

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

function RvalsGen(p::NamedTuple,ei::Electrode)
	N = ei.n*p.nsite
	R = Vector{Vector{Float64}}(undef,N)
	ix = -1
	if(ei.connectfrom=="-x")
		ix = -1
	elseif(ei.connectfrom=="+x")
		ix = p.nx+1
	else
		return 0 #something is busted
	end
	# offset to fix the 
	iRoffset = Int(0 + 0 + ix*p.nsite + ei.yrange[1]*p.nx*p.nsite + ei.zrange[1]*p.ny*p.nx*p.nsite)
        for iy = ei.yrange[1]:(ei.yrange[2]-1)
            for iz = ei.zrange[1]:(ei.zrange[2]-1)
                for isite = 0:(p.nsite-1)
                    iR = Int(1 + isite + p.nsite*(ix-ix + p.nx*((iy-ei.yrange[1]) + (iz-ei.zrange[1])*p.ny)))
                    #iR = Int(1 - iRoffset + isite + ix*p.nsite + iy*p.nx*p.nsite + iz*p.ny*p.nx*p.nsite)
                    #println("$ix $iy $iz $isite iR $iR")
                    ivec = Int.([ix,iy,iz,isite])
                    Rval = xyztor(p,ivec)
                    R[iR] = deepcopy(Rval)
                end
            end
	end
	return R # needs ⊗I(p.norb)⊗I(2) for full (spinful) hilbert space
end

# Defines a cᵦ†cₐ term 

function zeeman(Bvals::Vector{Vector{Float64}},  p::NamedTuple, ElectrodeInfo::Electrode)
	# only defined for S-like orbitals with lz = 0
	N = ElectrodeInfo.n*p.nsite*p.norb*2
	#N = p.n*p.nsite*p.norb*2
	#for i = 1:N
	zeeman = spzeros(ComplexF64, N, N)
	C = ħ/(2*m₀) #sans q factor -> eV
	#Bxvals = Bvals[:][1]
	for ax = 1:3
		BiVals = [B[ax] for B in Bvals]
		zeeman .+= 2*C*Diagonal(BiVals)⊗I(p.norb)⊗σ[ax]
	end
	return sparse(zeeman)
end


# return a vector of Σ(k) functions which return Σₖ(E) which return a sparse nsite × nsite matrix at a given energy
function genΣₖs(p::NamedTuple,ElectrodeInfo::Vector{Electrode})   
	nE = size(ElectrodeInfo)[1]
	Σks = Vector{Function}(undef,nE)
	for i = 1:nE # define a self-energy function for every electrode attached to device
		NNs = genNNs(p,ElectrodeInfo[i])
		P = changeBasis(p,ElectrodeInfo[i])
		Hslab, Hₗ, Hᵣ = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
		Σk(k) = Σgen(p,Hslab(k),Hₗ(k),Hᵣ(k),ElectrodeInfo[i],P)
		Σks[i] = Σk
	end
	return Σks
end

function xyzElectrodeSiteToI(ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	   ix = 0*ivec[1]; iy = ivec[2]; iz = ivec[3]
	return  ix + nx*(iy + ny*(iz)) + 1
end

function xyzElectrodeToI(p::NamedTuple, ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; iorb = ivec[4]
	return  iorb + p.norb*(ix + nx*(iy + ny*(iz))) + 1
end

function electrodeSiteToDeviceIndex(p::NamedTuple, ElectrodeInfo::Electrode,ivecContact::Vector{Int})
	# now construct an ivec for the site in the device
	if(ElectrodeInfo.connectfrom=="-x")
		ix = 0 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = (p.nx-1)
	end
	iy = ivecContact[2]+ElectrodeInfo.yrange[1]
	iz = ivecContact[3]+ElectrodeInfo.zrange[1]
	return ix + p.nx*(iy+p.ny*iz) + 1
end

function changeBasis(p::NamedTuple,ElectrodeInfo::Electrode)
	nE = ElectrodeInfo.n*p.nsite
	nD = p.n*p.nsite
	#println("nE = $nE; nD = $nD")
	Psite = spzeros(nD,nE)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	# only consider sites that will actually touch the slab
	ix = 0
	if(ElectrodeInfo.connectfrom=="-x")
		ix = p.nx-1 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = 0
	end
	for iy = 0:(ny-1)
		for iz = 0:(nz-1)
			ivec = [ix,iy,iz,0]
			contactSiteIndex = xyzElectrodeSiteToI(ElectrodeInfo,ivec)
			deviceSiteIndex = electrodeSiteToDeviceIndex(p,ElectrodeInfo,ivec)
			#println("Device site: $deviceSiteIndex, Contact site: $contactSiteIndex")
			Psite[deviceSiteIndex,contactSiteIndex] = 1
		end
	end
	return Psite⊗I(p.norb)⊗I(2)
end
			


function makeElectrodeH(p::NamedTuple,ElectrodeInfo::Electrode,edge_NNs::Vector{Hopping})
	kfilter = [0;1;1] # to ensure that self-energy in x is Γ centered
	function H(k::Vector{Float64})
		# hamiltonian describing the edges
		Hₑ = spzeros(ComplexF64, 2*p.nsite*p.norb*ElectrodeInfo.n, 2*p.nsite*p.norb*ElectrodeInfo.n)
		for NN in edge_NNs
			Δϕ = exp(im*(kfilter.*k)⋅(p.A*NN.N))
			#Δϕ = exp(im*k⋅(p.SLa₁*NN.N[1] + p.SLa₂*NN.N[2] + p.SLa₃*NN.N[3]))
			Hₑ[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]*Δϕ
			Hₑ[2*NN.b  , 2*NN.a-1] += NN.t[2,1]*Δϕ
			Hₑ[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]*Δϕ
			Hₑ[2*NN.b  , 2*NN.a  ] += NN.t[2,2]*Δϕ
		end
		return Hₑ
	end
	return H
end

function electrodeParams(p::NamedTuple,ElectrodeInfo::Electrode)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ep = (nx = nx, ny = ny, nz = nz, n = nz*nx*ny)
	return merge(p,ep)
end

function HcontactGen(p::NamedTuple,NNs::Vector{Hopping},ElectrodeInfo::Electrode)
	edgeNNs = Hopping[]
	LedgeNNs = Hopping[]
	RedgeNNs = Hopping[]
	nextLayerNNs = Hopping[]
	N = ElectrodeInfo.n*p.nsite*p.norb*2
	Rvals = RvalsGen(p,ElectrodeInfo)
	Bfield = ElectrodeInfo.A.(Rvals)
        Bfield = ElectrodeInfo.A.(Rvals)
	Hᵦ = zeeman(Bfield,p,ElectrodeInfo)
        #display(Hᵦ)
        H₀ = spzeros(ComplexF64,N, N) .+ Hᵦ
	#NNs = genNNs(p,ElectrodeInfo)
	pruneHoppings(NNs,p.prune) # cuts off the relevant matrix elements to make thin film
	for NN in NNs
		if(NN.N⋅[1;0;0] > 0)
			push!(RedgeNNs,deepcopy(NN))
		elseif(NN.N⋅[1;0;0] < 0)
			push!(LedgeNNs,deepcopy(NN))
		elseif(NN.edge==true)
			push!(edgeNNs,deepcopy(NN))
		else
			H₀[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H₀[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H₀[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H₀[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	Hc = makeElectrodeH(p,ElectrodeInfo,edgeNNs)
	Hₗ = makeElectrodeH(p,ElectrodeInfo,LedgeNNs)
	Hᵣ = makeElectrodeH(p,ElectrodeInfo,RedgeNNs)
	function Hslab(k::Vector{Float64})
		return Hc(k).+H₀	
	end	
	return Hslab, Hₗ, Hᵣ
end


#function Σgen(p::NamedTuple,H::Matrix,Hₗ::Matrix, Hᵣ::Matrix, ElectrodeInfo::Electrode, cutoff::Float64=10^-7*eV)
function Σgen(p::NamedTuple,H::SparseMatrixCSC,Hₗ::SparseMatrixCSC, Hᵣ::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-7*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    # so H needs to be instantiated and called outside of the loop
	H_coupling = Hₗ .+ Hᵣ# couples a layer to the infinite on both sides
    function Σ(E::Float64)
        Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*H_coupling'
        # converge the self energy 
        error = 1
        #for i = 1:40
        while error > cutoff
            #println("Σ convergence error loop: $error")
            Σ_guess0 = deepcopy(Σ_guess)
            Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- Σ_guess0)*H_coupling'
            error =  norm(Σ_guess.-Σ_guess0)
        end
        if(ElectrodeInfo.connectfrom=="-x")
                Σ_surf = Hᵣ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hᵣ'
                return P*Σ_surf*P'
        else
                if(ElectrodeInfo.connectfrom != "+x")
                        println("Something is very broken, check that your electrodes are in ±x")
                end
                Σ_surf = Hₗ*grInv((E+im*p.η)*I(n) .- H .- Σ_guess)*Hₗ'
                #show(size(Σ_surf))
                #show(size(P))
                return P*Σ_surf*P'
        end
    end
    return Σ
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
			H[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	return H, edgeNNs
end

# Add a bond to the list of bonds, given some list of bonds, coefficient in spin basis, index of both sites, and param list
function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0],""))
	push!(NNs,NN)
end

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


function weylHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        t = 3*nextsite(iorb)*p.t*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; di[5] = nextsite(iorb); ib = Int.(ia + di)
			Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
			t = (-im/2)*dir*p.t*σ[ax]
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
	end
end
	
function nextsite(isite::Int)
	return -2*isite + 1
end

function genNNs(p,ElectrodeInfo::Electrode) # all of the terms in the hamiltonian get added here, get back the relevant bonds
	n = p.n
	NNs = Hopping[]
	ep = electrodeParams(p,ElectrodeInfo) # electrodeparams
	ix = 0
	nx = 1
	for iy = 0:(ep.ny-1)
		for iz = 0:(ep.nz-1)
		     for isite = 0:(ep.nsite-1)
			for iorb = 0:(ep.norb-1)
                            ia = (copy(Int.([ix,iy,iz,isite,iorb])));
                            if(ElectrodeInfo.type=="weyl")
                                    weylHopping(ep,NNs,ia)
                            end
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
		pib = ib - [mod(ib[1],nx),mod(ib[2],ep.ny),mod(ib[3],ep.nz)]
		#pib = ib - [mod(ib[1],p.nx),mod(ib[2],p.ny),mod(ib[3],p.nz)]
		if(pib⋅pib != 0) # if vector is distinctly outside of 1st lattice
			NN.N = Int.([round(pib[1]/(nx)),round(pib[2]/ep.ny),round(pib[3]/ep.nz)])
			NN.b = xyztoi(ep,NN.ib, NN.N)
			NN.edge = true
			#println("$(NN.N)")
		end
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
