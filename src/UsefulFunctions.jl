
module UsefulFunctions

export rotate, average, ⊗, recInv

using LinearAlgebra
using SparseArrays

⊗(A,B) = kron(A,B)
×(u,v) = cross(u,v)
#CartProd(x,y) = [[x,y] for i in x for j in y]
function CartProd(Vecs)
        type = typeof(Vecs[1][1])
        #init = Any[]
        #nextstep = Any[]
        init = copy(Vecs[1])
        #nextstep = deepcopy(init)
        for iarg in 2:size(Vecs)[1]
            #v1 = Vecs[iarg-1]
            v2 = Vecs[iarg]
            nextstep = [vcat(i,j) for i in init for j in v2]
            init = copy(nextstep)
        end
        return init
end

function mkfolder(path)
	if(isdir(path))
		println("$path already exists...")
		#rm(path, recursive=true)
	else
		mkdir(path)
	end
end

function b(row::Int,B::Int=4)
    return ((row-1)*B+1):(row*B)
end

function buildM!(Mrows::Vector{Int},Mcols::Vector{Int},elems::Vector{ComplexF64}, rows::UnitRange{Int},cols::UnitRange{Int},submat)
    offsetrow = rows[1]; offsetcol = cols[1]
    for row in rows
        for col in cols
            push!(Mrows,row); push!(Mcols,col)
            push!(elems,submat[row-offsetrow+1,col-offsetcol+1])
        end
    end
end    

function pGrInv(M::SparseMatrixCSC, B::Int=4, offdiag = true)
    n = size(M)[1];
    nBs = Int(n/B) # number of diagonal blocks
    LCblocks = Matrix[] # go down the diagonal from left
    RCblocks = Matrix[] # go down the diagonal from right 
    # see Gerhardt Klimeck's powerpoint
    Minv = spzeros(ComplexF64,n,n)
    rows = Int[]; cols = Int[]; elems = ComplexF64[];
    gʳL₀ = inv(Array(M[1:B,1:B]))
    push!(LCblocks,gʳL₀)
    for ib = 2:nBs
        gʳLprior = LCblocks[ib-1]
        # off-diagonal block coupling hamiltonian
        V = M[b(ib-1,B),b(ib,B)]
        #display(V)
        # on-diagonal block hamiltonian
        Dᵢ = M[b(ib,B),b(ib,B)]
        # effective surface green's function including this point
        gʳLᵢ = inv(Array(Dᵢ - V'*gʳLprior*V))
        push!(LCblocks,copy(gʳLᵢ))
    end
    # now we will go back up the diagonal and incorporate coupling from right
    Gʳend = last(LCblocks)
    Gʳᵢ = Gʳend
    buildM!(rows,cols,elems,b(nBs,B),b(nBs,B),Gʳend)
    Gʳplus = copy(Gʳend)
    if offdiag == "legacy"
        for ib = reverse(1:(nBs-1))
            gʳLᵢ = LCblocks[ib]
            # off-diagonal block coupling hamiltonian
            V = M[b(ib,B),b(ib+1,B)]
            # see Klimeck's pwpt
            Gʳᵢ = gʳLᵢ*(I(B) + V*Gʳplus*V'*gʳLᵢ)
            # add one extra diagonal
            Gʳᵢoffdiag = Array(-gʳLᵢ*V*Gʳplus)
            buildM!(rows,cols,elems,b(ib+1,B),b(ib,B),transpose(Gʳᵢoffdiag))
            buildM!(rows,cols,elems,b(ib,B),b(ib+1,B),Gʳᵢoffdiag)
            # now do the block-diagonal Gʳ
            buildM!(rows,cols,elems,b(ib,B),b(ib,B),Array(Gʳᵢ))
            Gʳplus = copy(Gʳᵢ)
        end
    elseif offdiag == "transport"
        # go back up the diagonal and also do bottom row
        Gʳbotplus = copy(Gʳend)
        for ib = reverse(1:(nBs-1))
            gʳLᵢ = LCblocks[ib]
            # off-diagonal block coupling hamiltonian
            V = M[b(ib,B),b(ib+1,B)]
            # see Klimeck's pwpt
            # now do the block-diagonal Gʳ
            Gʳᵢ = gʳLᵢ*(I(B) + V*Gʳplus*V'*gʳLᵢ)
            buildM!(rows,cols,elems,b(ib,B),b(ib,B),Array(Gʳᵢ))
            Gʳᵢbot = transpose(Array(transpose(-LCblocks[ib])*V*transpose(Gʳbotplus)))
            buildM!(rows,cols,elems,b(ib,B),b(nBs,B),transpose(Gʳᵢbot))
            buildM!(rows,cols,elems,b(nBs,B),b(ib,B),Gʳᵢbot)
            Gʳplus = copy(Gʳᵢ)
            Gʳbotplus = copy(Gʳᵢbot)
        end
        # generate the right-connected gʳs
        gʳR₀ = inv(Array(M[b(nBs,B),b(nBs,B)]))
        push!(RCblocks,copy(gʳR₀))
        gʳRprior = gʳR₀
        for ib = reverse(1:(nBs-1))
            # off-diagonal block coupling hamiltonian
            V = M[b(ib,B),b(ib+1,B)]
            #display(V)+
            # on-diagonal block hamiltonian
            Dᵢ = M[b(ib,B),b(ib,B)]
            # effective surface green's function including this point
            gʳRᵢ = inv(Array(Dᵢ - V*gʳRprior*V'))
            push!(RCblocks,copy(gʳRᵢ))
            gʳRprior = gʳRᵢ
        end
        # and now do the top row
        reverse!(RCblocks)
        Gʳtopmin = Gʳᵢ
        for ib = 2:nBs
            gʳRᵢ = RCblocks[ib]
            V = M[b(ib,B),b(ib-1,B)]
            Gʳᵢtop = -transpose(gʳRᵢ*V*transpose(Gʳtopmin))
            buildM!(rows,cols,elems,b(1,B),b(ib,B),Gʳᵢtop)
            buildM!(rows,cols,elems,b(ib,B),b(1,B),transpose(Gʳᵢtop))
            Gʳtopmin = copy(Gʳᵢtop)
        end
    else # or for case without the off-diagonal blocks
        for ib = reverse(1:(nBs-1))
            gʳLᵢ = LCblocks[ib]
            # off-diagonal block coupling hamiltonian
            V = M[b(ib,B),b(ib+1,B)]
            # see Klimeck's pwpt
            Gʳᵢ = gʳLᵢ*(I(B) + V*Gʳplus*V'*gʳLᵢ)
            buildM!(rows,cols,elems,b(ib,B),b(ib,B),Array(Gʳᵢ))
            Gʳplus = copy(Gʳᵢ)
        end
    end
    return sparse(rows,cols,elems)
end

function int(i::Float64)
    return Int(round(i,sigdigits=1))
end

#CartProd(x,y) = [repeat(x, inner=[size(y,1)]) repeat(y, outer=[size(x,1)])]
#CartProd(x,y,z) = [[x,y,z] for i in x, for j in y, for ]
function average(v)
	N = size(v)[1]
	return sum(v)*(1/N)
end

avg(v) = average(v)

function genTetBZ(p::NamedTuple,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    X1 = p.kdict["X₁"];
    X2 = p.kdict["X₂"];
    X3 = p.kdict["X₃"];
    function divFixNaN(a::Int,b::Int) # for this particular instance, n/0 represents a Γ-centred sampling @ k = 0. 
            if(b==0)
                    return 0
            else
                    return a/b
            end
    end
    for ix = -nx:nx
        for iy = -ny:ny
            for iz = -nz:nz
                kindex = [iy + ny + 1; iz + nz + 1]
                k = divFixNaN(ix,nx)*X1 + divFixNaN(iy,ny)*X2 + divFixNaN(iz,nz)*X3
                kweight = 1
                if(abs(ix) == nx)
                    kweight *= 1/2
                end
                if(abs(iy) == ny)
                    kweight *= 1/2
                end
                if(abs(iz) == nz)
                    kweight *= 1/2
                end
                push!(kpoints,k)
                push!(kindices,kindex)
                push!(kweights,kweight)
            end
        end
    end
    ksum = sum([w for w in kweights])
    kweights = (1/ksum).*kweights
    return kpoints, kweights, kindices
end

#grInv(A) = inv(Array(A))
exactInv(A) = inv(Array(A))

# using Woodbury matrix recursion identity
function recInv(M::Union{SparseMatrixCSC, SubArray}, sizecutoff::Int=64)
    nx = size(M)[2]
    ny = size(M)[1]
    if(nx <= sizecutoff || ny <= sizecutoff)
        return sparse(exactInv(M))
    else
        hxu = Int(round(nx/2))+1; hyu = Int(round(ny/2))+1
        hxl = hxu -1; hyl = hyu -1
        #A = deepcopy(M[1:hyl,1:hxl]); U = deepcopy(M[1:hyl,hxu:nx]);
        #V = deepcopy(M[hyu:ny,1:hxl]); C = deepcopy(M[hyu:ny,hxu:nx]);
        A = M[1:hyl,1:hxl];  U = M[1:hyl,hxu:nx];
        V = M[hyu:ny,1:hxl]; C = M[hyu:ny,hxu:nx];
        Cinv = recInv(C,sizecutoff)
        Γ = Cinv*V; Δ = U*Cinv
        #display(typeof(A))
        #display(typeof(sparse(deepcopy(Δ*V))))
        #display(typeof(sparse(A)-Δ*sparse(V)))
        #display(typeof(A.-deepcopy(Δ*V)))
        Σ = recInv(A .- Δ*V,sizecutoff)
        Minv = [Σ          -Σ*Δ;
                -Γ*Σ   Cinv .+ Γ*Σ*Δ]
        
        #=
        Cinv = sE(recInv(C,sizecutoff))
        Γ = Cinv*sE(V); Δ = sE(U)*Cinv
        Σ = sE(recInv(A .- Δ*V,sizecutoff))
        Minv = [Σ          -Σ*Δ;
                -Γ*Σ  (Cinv .+ Γ*Σ*Δ)]
        =#
        #Minv[1:hyl,1:hxl] .= Σ
        #Minv[1:hyl,hxu:nx] .= -Σ*U*Cinv
        #Minv[hyu:ny,1:hxl] .= -Cinv*V*Σ
        #Minv[hyu:ny,hxu:nx] .= Cinv .+ Cinv*V*Σ*U*Cinv
        #if(typeof(Minv) == SparseMatrixCSC)
        #    return dropzeros(Minv)
        #else
        return Minv
        #end
    end
end



grInv = recInv

function sparseEnough(M::Union{SparseMatrixCSC,SubArray})
       return M
       n = size(M)[1]*size(M)[2]
       nz = nnz(M)
       sparsitycutoff = nz
       if(nz/n^2 < sparsitycutoff)
           return dropzeros(M)
       else
           return Array(M)
       end
end

sE = sparseEnough

findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]


mutable struct Atom
	name
	R
	orbitals # vector with -1,0,1 for low/relevant/high
	atomNum
	valenceNum
end

mutable struct Connection
	Atom1
	Atom2 
	NN # number denoting collection of vectors to use in lattice
	radiusCoeff # scaling coefficient to put on vectors
#end

function ∇(f::Function, R, dx=10^-11)
	ndim = size(R)[1]
	grad = zeros(ndim)
	for i = 1:ndim
		dR = zeros(ndim)
		dR[i] = dx
		grad[i] = (f(R+dR)-f(R-dR))/(2*dx)
	end
	return grad
end

rot(θ) = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]


function Hessian(f::Function, R, dx = 10^-11)
	ndim = size(R)[1]
	hess = zeros(ndim,ndim)
	for i = 1:ndim
		for j = 1:ndim
			dR₁ = zeros(ndim)
			dR₁[i] = dx
			dR₂ = zeros(ndim)
			dR₂[j] = dx
			hess[i,j] = (f(R + dR₁ + dR₂) - f(R + dR₁ - dR₂) - f(R - dR₁ + dR₂) + f(R - dR₁ - dR₂))/(4*dx^2)
		end
	end
	return hess
end

function rotate(θ::Float64)
	return [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
end




function MtoA(M)
	return [M[i,:] for i in 1:size(M)[1]]
end
#=function eig(A, conv)
	nA = size(A)[1]
	n = nA # num of eigvls to calculate
	λ = zeros(Float64,n)
	d = 100 # value added to diagonal to make all eigvls positive
	vecs = zeros(ComplexF64,n,n)
	M = zeros(ComplexF64,nA,nA)
	M .= A + d*I(nA)
	for i in 1:n
		#get initial random vectors
		v = rand(ComplexF64,n)
		v .= v/norm(v)
		v_prior = rand(ComplexF64,n)

		#will find the largest eigenvector/eigvl of M
		error = 100
		#println("eigvl # $i")
		#println("\n Eigvl #$i")
		#j = 0
		#for Ω = 1:10000
		while(error > conv)
			v .= M*v
			#println("error = $error")
			error = abs(1 - norm(v)/norm(v_prior)) #works off eigenvalue convergence, instead of eigvec convergence
			v_prior .= v
			#display(v)
			v .= v/norm(v)
			#j += 1
		end
		E = norm(M*v)

		#subtract off contribution of Eₙ
		M .-=  E*v*v'

		#write into output 
		λ[i] = E - d
		for j in 1:n
			vecs[i,j] = v[j]
		end
		#println("norm(A - PDP†) = \n")
		#display(norm(A .-vecs'*Diagonal(λ)*vecs))
		#println("")
	end
	return λ, vecs
end =#

for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
       end

end

end
