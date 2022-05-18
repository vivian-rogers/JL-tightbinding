
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
S² = [ 
      (1/2)*(1/2 + 1) 	0
      0			(1/2)*(1/2 + 1)
     ]


σ = Vector{Matrix}(undef,3)
σ[1] = σ₁; σ[2] = σ₂; σ[3] = σ₃

S₁ = (1/2)*σ₁; S₂= (1/2)*σ₂; S₃ = (1/2)*σ₃; 
S = cat(S₁,S₂,S₃, dims = 3);



#	s	py	px	pz	dx2-y2	dyz	dxy	dz2	dzx
L2 =	[0, 	1,	1,	1,	2,	2,	2,	2,	2]
Lz =	[0, 	1, 	0, 	-1,	2,	1,	0,	-1,	-2]

L² = Diagonal([L2[i]*(L2[i]+1) for i = 1:9])
L₃ = Diagonal(Lz)

Lplus = [√( L2[i]*(L2[i]+1) - Lz[i]*(Lz[i]+1)) for i = 2:9]
Lminus = [√( L2[i]*(L2[i]+1) - Lz[i]*(Lz[i]-1)) for i = 1:8]

L₊ = Tridiagonal(zeros(8),zeros(9),Lplus)
L₋ = Tridiagonal(Lminus,zeros(9),zeros(8))


L₁ = Tridiagonal((1/2)*(L₊ .+ L₋))

L₂ = Tridiagonal((1/2)*(-im*L₊ .+ im*L₋)) 

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
