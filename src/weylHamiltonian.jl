

#push!(LOAD_PATH, "./src/")
module nearestNeighbors

using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using SparseArrays
using Operators
using GMaterialParameters
using Bands



function itot(ix,iy,iz) 
	# periodic boundary conditions
	ix = mod(ix-1,nx)+1; iy = mod(iy-1,ny)+1; iz = mod(iz-1,nz)+1
	return ix + (nx-1)*iy + (nx-1)*(ny-1)*iz
end


function Rvals(nᵢ, aᵢ)
	# generates the list of R values for each site
	nx = nᵢ[1]; ny = nᵢ[2]; nz = nᵢ[3]
	N = nx*ny*nz
	#println("N = $N\n")
	#δ = (a/sqrt(3))*[cosd(30);sind(30)];
	Rvals = zeros(N,2)
	function R(ix,iy,iz)
		return (ix-1)*aᵢ[1]*[1;0;0] + (iy-1)*aᵢ[2]*[0;1;0] + (iz-1)*aᵢ[3]*[0;0;1]
	end
	for ix = 1:nx
		for iy = 1:ny
			for iz = 1:nz
				Rvals[itot(ix,iy,iz),:] = R(ix,iy,iz)
			end
		end
	end
	return Rvals
end

function Rlist(nᵢ,Rarr) # puts list in different format if needed
	N = nᵢ[1]*nᵢ[2]*nᵢ[3]
	#N = 2*λ^2
	RvalList = Rarr
	#RvalList = Rvals(N)
	#Rarr = Array[]
	#for i = 1:N
	#	Rarr = Rs[i,:]
	#end
	Rs = [RvalList[i,:] for i in size(RvalList)[1]]
	return Rs 
end

function Rofg(g,λ) 
	# grid-pos to position
	x = g[1]; y = g[2];
	δ = (a/sqrt(3))*[cosd(30);sind(30)];
	a₁ = a*[1;0]; a₂ = a*[cosd(60);sind(60)]; 
	R = floor((x-1)/2)*a₁+(y-1)*a₂
	if(mod(x,2) == 0)
		R += δ
	end
	return R
end


function g2i(x,y,λ) # grid-pos to index in hamiltonian
	return x + (y-1)*(2*λ)
	# i = x + (y-1)*2λ; (i - x)/2λ + 1= y 
	# i - (y-1)*2λ = x
end


function i2g(i,λ) # reverse above
	y = ceil(i/(2*λ))
	x = i - (y-1)*(2λ)
	return x,y
end



# defined a bond, with the site-indices, radii, etc, hopping value. Can be modified. 
mutable struct Hopping
	a::Int #site 1
	b::Int #site 2
	r # radius from a to b
	gA # grid-value of point a
	gB # grid-value of point B
	t  # hopping parameter affiliated with c†₂c₁
	edge::Bool
	N # vector describing the [n₁;n₂]⋅[a₁;a₂] superlattice unit cell
end

# turns the hopping structs into a matrix
function nnHoppingMat(NNs,nᵢ)
	N = nᵢ[1]*nᵢ[2]*nᵢ[3]
	H = spzeros(ComplexF64,N,N)⊗I(2)
	edgeNNs = Any[]
	for NN in NNs
		if(NN.edge == true)
			push!(edgeNNs,deepcopy(NN))
		else
			Hi = spzeros(N,N)
			Hi[NN.a,NN.b] = 1
			H .+= Hi[NN.a,NN.b]⊗NN.t
			#H[NN.a,NN.b] = NN.t
		end
	end
	return H, edgeNNs
end

function nn1(nᵢ,Rvals,t₀)
	# generates the nearest neighbors for the weyl hamiltonian	
	nx = n
	N = nᵢ[1]*nᵢ[2]*nᵢ[3]
	NNs = Any[]
	for ix = 1:(nx)
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
	
	#=for y = 1:λ
		for x = 1:2:(2*λ) #site A
			i1 = g2i(x,y,λ);
			r1 = [x;y]
		
			r2 = [x+1,y-1] #down one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)

			r2 = [x-1,y] #left one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x+1,y] #right one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)
		end
		for x = 2:2:(2*λ) #site B
			i1 = g2i(x,y,λ);
			r1 = [x;y]
			
			r2 = [x-1,y+1] #up one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x-1,y] #left one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x+1,y] #down one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,t₀,false,[0;0])
			push!(NNs,NN)
		end
	end
	for NN in NNs # let us fix the edge designation
		x = NN.gB[1]; y = NN.gB[2] # gives grid-value of the site off the edge
		#show(NN)
		if(x<1 || x>(2*λ) || y<1 || y>λ)
			
			# generate basis with wraparound
			xnew = mod(x-1,2*λ)+1
			ynew = mod(y-1,λ)+1
			NN.b =g2i(xnew,ynew,λ)
			
			# fill in N vector
			
			if(y<1)
				ny = -1
			elseif(y>λ)
				ny = 1
			else
				ny = 0
			end
			
			if(x<1)
				nx = -1
			elseif(x>2*λ)
				nx = 1
			else
				nx = 0
			end
			NN.N = [nx;ny]
			NN.edge = true
		end
	end =#
	return NNs
end


function nn2(λ,Rvals,U₂)
	
	# this is somewhat misleading -- these "hoppings" are not actually use for hopping. Hijack the hopping struct for the U₂ NNN repulsion
	N = 2*λ^2	
	NNs = Any[]
	r2_dg = [[0,-1],[2,-1],[0,2],[0,1],[-2,1],[-2,0]] # r2 grid offsets
	for y = 1:λ
		for x = 1:(2*λ) #site A
			i1 = g2i(x,y,λ);
			r1 = [x;y]
			#do NNs
			for i = 1:6
				dx =r2_dg[i][1]; dy = r2_dg[i][2]
				r2 = [x+dx,y+dy]
				i2 = g2i(r2[1],r2[2],λ);
				NN = Hopping(i1,i2,radius(r1,r2,λ),r1,r2,U₂,false,[0;0])
				push!(NNs,NN)
			end
		end
	end
	for NN in NNs # let us fix the edge designation
		x = NN.gB[1]; y = NN.gB[2] # gives grid-value of the adjoint
		#show(NN)
		if(x<1 || x>(2*λ) || y<1 || y>λ)
			
			# generate basis with wraparound
			xnew = mod(x-1,2*λ)+1
			ynew = mod(y-1,λ)+1
			NN.b =g2i(xnew,ynew,λ)
			
		end
	end
	return NNs
end


function radius(g1,g2,λ) # distance vector between two grid-sites
	r1 = Rofg(g1,λ)
	r2 = Rofg(g2,λ) 
	return (r2 - r1)
end

for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end

