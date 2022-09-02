
module Operators

using LinearAlgebra
using Constants
#using SparseArrays

⊗(A,B) = kron(A,B)
×(u,v) = cross(u,v)

τ₀ = [
      1 0
      0 1
     ]

τ₁ = [
      0  1 
      1  0
     ] 

τ₂ = [
      0 -im 
      im  0
     ] 

τ₃ = [
      1  0
      0 -1
     ]

σ₀ = [
      1 0
      0 1
     ]

σ₁ = [
      0  1 
      1  0
     ] 
σ₂ = [
      0 -im 
      im  0
     ] 
σ₃ = [
      1  0
      0 -1
     ]


S1z =             [1 0 0;
                   0 0 0;
                   0 0 -1]

S1y =    sqrt(1/2)*[0 -im 0;
                    im 0 -im;
                    0 im 0]

S1x =   sqrt(1/2)*[0 1 0;
                   1 0 1;
                   0 1 0]


S² = [ 
      (1/2)*(1/2 + 1) 	0
      0			(1/2)*(1/2 + 1)
     ]


σ = Vector{Matrix}(undef,3)
σ[1] = σ₁; σ[2] = σ₂; σ[3] = σ₃



S₁ = (1/2)*σ₁; S₂= (1/2)*σ₂; S₃ = (1/2)*σ₃; 
S = cat(S₁,S₂,S₃, dims = 3);

function matdot(A::Vector{Matrix},B::Vector{Float64})
    sum = zeros(ComplexF64,size(B[1]))
    for i in eachindex(B)
        sum .+= B[i]*A[i]
    end
    return sum
end

function matdot(A::Vector,B::Vector{Matrix})
    sum = zeros(ComplexF64,size(B[1]))
    for i in eachindex(A)
        sum .+= A[i]*B[i]
    end
    return sum
end

function distFromDW(p,Rvals::Vector{Vector{Float64}})
        a = p.A[1,1]
        DWs = a*[0,0.5,1] # x pos of domain walls
        nearestDW = [minimum([abs(R[1] - DWpos) for DWpos in DWs]) for R in Rvals]

        return Diagonal(1/nm*nearestDW)
end

function zpos(Rvals::Vector{Vector{Float64}})
	zmax = maximum([R[3] for R in Rvals])
	return Diagonal([(R[3]-zmax)/nm for R in Rvals])
end



# see: Spherical change of basis, complex basis wikipedia
pmToOrb = sqrt(1/2)*[-1 0 1; im 0 im; 0 sqrt(2) 0]
#pmToOrb = sqrt(1/2)*[1 0 1; im 0 -im; 0 sqrt(2) 0]

L₁ = pmToOrb*S1x*pmToOrb'
L₂ = pmToOrb*S1y*pmToOrb'
L₃ = pmToOrb*S1z*pmToOrb'

L = Vector{Matrix}(undef,3)

L[1] = L₁; L[2] = L₂; L[3] = L₃

Ldotσ = sum([L[i]⊗σ[i] for i = 1:3])

#J₁ = L₁⊗I(2) .+ I(9)⊗S₁
#J₂ = L₂⊗I(2) .+ I(9)⊗S₂
#J₃ = L₃⊗I(2) .+ I(9)⊗S₃

#J² = J₁*J₁ .+ J₂*J₂ .+ J₃*J₃
#altLdotS = (1/2)*(J² .- (L²⊗I(2) +  I(9)⊗S²))

#LdotS = (L₁⊗I(2))*(I(9)⊗S₁) + (L₂⊗I(2))*(I(9)⊗S₂) + (L₃⊗I(2))*(I(9)⊗S₃)


function gradients(nx,ny,nz, aᵢ, pbx = 0, pby = 0, pbz = 0) # generates gradient matrices for square lattice
	∇₁ = ∇₂ = ∇₃ = zeros(nx*ny*nz)
	∇ᵢ = [∇₁; ∇₂; ∇₃]
	function itot(ix,iy,iz) 
		# periodic boundary conditions
		ix = mod(ix-1,nx)+1; iy = mod(iy-1,ny)+1; iz = mod(iz-1,nz)+1
		return ix + (nx-1)*iy + (nx-1)*(ny-1)*iz
	end
	# instantiate gradient at all points in lattice 
	for ix = 1:(nx+pbx-1)
		for iz = 1:(ny+pby-1)
			for iz = 1:(nz+pby-1)
				i = itot(ix,iy,iz)
				# generate grad x,y,z
				for dir = eachindex(∇ᵢ)
					# central derivative
					for bond = [-1,1]
						ivec = [ix,iy,iz]
						ivec[dir] += bond
						ibond = itot(ivec[1], ivec[2], ivec[3])
						∇ᵢ[dir][i,ibond] = bond/(2*aᵢ[dir])
					end
				end
			end
		end
	end
	return ∇ᵢ[1], ∇ᵢ[2], ∇ᵢ[3]
end

function gradients(nx,ny,nz, aᵢ, pbx = 0, pby = 0, pbz = 0) # generates gradient matrices for square lattice
	∇₁ = ∇₂ = ∇₃ = zeros(nx*ny*nz)
	∇ᵢ = [∇₁; ∇₂; ∇₃]
	function itot(ix,iy,iz) 
		# periodic boundary conditions
		ix = mod(ix-1,nx)+1; iy = mod(iy-1,ny)+1; iz = mod(iz-1,nz)+1
		return ix + (nx-1)*iy + (nx-1)*(ny-1)*iz
	end
	# instantiate gradient at all points in lattice 
	for ix = 1:(nx+pbx-1)
		for iz = 1:(ny+pby-1)
			for iz = 1:(nz+pby-1)
				i = itot(ix,iy,iz)
				# generate grad x,y,z
				for dir = eachindex(∇ᵢ)
					# central derivative
					for bond = [-1,1]
						ivec = [ix,iy,iz]
						ivec[dir] += bond
						ibond = itot(ivec[1], ivec[2], ivec[3])
						∇ᵢ[dir][i,ibond] = bond/(2*aᵢ[dir])
					end
				end
			end
		end
	end
	return ∇ᵢ[1], ∇ᵢ[2], ∇ᵢ[3]
end

Hₙ = Diagonal([-Ry*q^2/n^2 for n = 1:6])


#function L₃(l) 
#	return Diagonal([ħ*m for m = -l:l])
#end


for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
