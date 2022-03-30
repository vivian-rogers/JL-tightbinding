module MaterialParameters

using Constants
using LinearAlgebra
#a = 1
# for ScN


struct Atom
	name # name of atom, ex: "Co1:
	type # the type of atom, ex: "Fe"
	i # index, for whatever system
	R # location of atom in metres
	#NNs # list of nearest neighbors from this atom
end

a = 6.2*Å
#δ = (1/√(3))*a
t = 2.8*eV
ε = 0*eV

a₁ = (a/2)*[0;1;1]; a₂ = (a/2)*[1;0;1]; a₃ = (a/2)*[1;1;0];
A = hcat(a₁, a₂, a₃);
#show(A)
B = transpose((2*π)*inv(A));

#reciprocal lattice
b₁ = B[:,1]; b₂ = B[:,2]; b₃ = B[:,3];

Atoms = Array{Atom, 1}(undef, 4)


kdict = Dict(
	    "Γ" => B*[0;	0;	0],
	    "K" => B*[3/8; 	3/8;	3/4],
	    "U" => B*[5/8; 	1/4;	5/8],
	    "L" => B*[1/2;	1/2;	1/2],
	    "W" => B*[1/2;	1/4;	3/4],
	    "X" => B*[1/2;	0;	1/2],
	    )	    


Atoms[1] = Atom("Co1","Co",1,0*(a₁ + a₂ + a₃))
Atoms[2] = Atom("Sc","Sc",2,0.25*(a₁ + a₂ + a₃))
Atoms[3] = Atom("Co2","Co",3,0.50*(a₁ + a₂ + a₃))
Atoms[4] = Atom("Sb","Sb",4,0.75*(a₁ + a₂ + a₃))


ϵ_onsite = eV*[-0.5,0.2,-0.25,-3]
tdict = Dict(
	     ["Co","Sc"] => 0.4*eV,
	     ["Sc","Co"] => 0.4*eV,
	     ["Co","Sb"] => 0.5*eV,
	     ["Sb","Co"] => 0.5*eV,
	     ["Sc","Sb"] => 1*eV,
	     ["Sb","Sc"] => 1*eV,
	     ["Co","Co"] => 0.2*eV,
)

#nn1_distance = a/2;
#nn2_distance 


#A_uc = a^2*sind(60) # area of the unit cell


#real-space lattice vectors
function kdictGen(λ,θ=0)
rot = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]

a₁ = rot*λ*a*[1;0]; a₂ = rot*λ*a*[cosd(60);sind(60)]; 
A = hcat(a₁, a₂);
B = transpose((2*π)*inv(A));

#reciprocal lattice
b₁ = B[:,1]; b₂ = B[:,2];

kdict = Dict(
	    "γ" => B*[0;     0],
	    "κ" => B*[2/3; 1/3],
	    "κ'" => B*[-2/3; -1/3],
	    "κ'2" => B*[1/3; -2/3],
	    "μ" => B*[1/2; 0],
	    "μ'" => B*[-1/2; 0],
	    )	    
return kdict
end

function rLattice(λ)
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	A = hcat(a₁, a₂);
	B = transpose((2*π)*inv(A));

	#reciprocal lattice
	b₁ = B[:,1]; b₂ = B[:,2];
	return b₁, b₂
end

function kgrid2(n,λ)
	b₁, b₂ = rLattice(λ)
	kdict = kdictGen(λ)
	klist = Any[]
	weights = Any[]
	#γ = kdict["γ"]
	#κp = kdict["κ'2"]
	#κ = kdict["κ"]
	if(n==0) #gamma centred
		push!(klist,[0,0])
		push!(weights,1)
	elseif(n==1) # do the corners of the IBZ
		push!(klist,0*b₁+0*b₂)
		push!(weights,2/3)
		push!(klist,0*b₁+1*b₂)
		push!(weights,1/3)
		#push!(klist,0*b₁+0*b₂)
		#push!(weights,1/4)
		#push!(klist,1*b₁+1*b₂)
		#push!(weights,1/3)
	else # evenly spaced grid over IBZ triangle
		# corners
		dk = 1/n
		#κp = kdict["κ'"] 
		#κ = kdict["κ"] 
		#γ = kdict["γ"] 
		g = 0*b₁ + 0*b₂
		g2 = 1*b₁ + 1*b₂
		push!(klist,g)
		push!(weights,1/3)
		push!(klist,b₂)
		push!(weights,1/6)
		push!(klist,b₁)
		push!(weights,1/6)
		push!(klist,g2)
		push!(weights,1/3)
		# do edges of IBZ
		for i = (dk):dk:(1-dk)
			k = (i)*g + (1-i)*b₁
			push!(klist,k)
			push!(weights,1/2)
			k = (i)*b₁ + (1-i)*g2
			push!(klist,k)
			push!(weights,1/2)
			k = (i)*g2 + (1-i)*b₂
			push!(klist,k)
			push!(weights,1/2)
			k = (i)*b₂ + (1-i)*g
			push!(klist,k)
			push!(weights,1/2)
		end	
		# do interior of IBZ
		# imagine tiling a ternary plot with [a,b,c]⋅[γ,κ,κp]
		for i = dk:dk:(1-dk)
			for j = dk:dk:(1-dk)
				k = i*b₁ + j*b₂
				push!(klist,k)
				push!(weights,1)
			end
		end

	end
	wsum = sum(weights)
	weightlist = [w/wsum for w in weights]
	return klist, weightlist
end

function kgrid(n,λ)
	kdict = kdictGen(λ)
	klist = Any[]
	weights = Any[]
	γ = kdict["γ"]
	κp = kdict["κ'2"]
	κ = kdict["κ"]
	if(n==0) #gamma centred
		push!(klist,γ)
		push!(weights,1)
	elseif(n==1) # do the corners of the IBZ
		push!(klist,γ)
		push!(weights,1/3)
		push!(klist,κp)
		push!(weights,1/3)
		push!(klist,κ)
		push!(weights,1/3)
	else # evenly spaced grid over IBZ triangle
		# corners
		dk = 1/n
		#κp = kdict["κ'"] 
		#κ = kdict["κ"] 
		#γ = kdict["γ"] 
		push!(klist,κp)
		push!(weights,1/6)
		push!(klist,κ)
		push!(weights,1/6)
		push!(klist,γ)
		push!(weights,1/6)
		# do edges of IBZ
		for i = (dk):dk:(1-dk)
			k = (i)*κp + (1-i)*γ
			push!(klist,k)
			push!(weights,1/2)
			k = (i)*κ + (1-i)*κp
			push!(klist,k)
			push!(weights,1/2)
			k = (i)*γ + (1-i)*κ
			push!(klist,k)
			push!(weights,1/2)
		end	
		# do interior of IBZ
		# imagine tiling a ternary plot with [a,b,c]⋅[γ,κ,κp]
		dt = 1/(n)
		for a = dt:dt:(1-2*dt)
			for b = dt:dt:(1-dt-a)
				c = 1 - a - b
				k = a*γ + b*κ + c*κp
				push!(klist,k)
				push!(weights,1)
			end
		end

	end
	wsum = sum(weights)
	weightlist = [w/wsum for w in weights]
	return klist, weightlist
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








#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
