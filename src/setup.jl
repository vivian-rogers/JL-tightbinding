
module UsefulFunctions

using LinearAlgebra

#export ⊗,⋅,×,Atom

⊗(A,B) = kron(A,B)
×(u,v) = cross(u,v)
CartProd(x,y) = [(x,y) for i in x for j in y]


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

function nn1_vects(a)
	vects = Any[]
	for i = 1:3; 
		for j = [-1,1]
			# r₁ vector is defined in src/HolmiumParameters, may warrant reformulation in terms of lattice vectors	
			r = zeros(3);
			r[i] = j*(a/2);
			push!(vects,r)
		end
	end
	return vects
end


function nn2_vects(a)
	vects = Any[]
	for i = 1:3; 
		for j = [-1,1]
			for k = [-1,1]
				r = ones(3);
				#just consider only the two relevant dimensions
				dims = [1,2,3]
				dims = filter!(dim->dim!=i,dims)
				r[i] = 0
				r[dims[1]] = j*a
				r[dims[2]] = k*a
				push!(vects,r)
			end
		end
	end
	return vects
end

rot(θ) = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]

function avg(v)
	N = size(v)[1]
	sum = 0
	for i in v
		sum += i
	end
	return sum/N
end

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
