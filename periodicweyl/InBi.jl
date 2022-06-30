push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module InBi
using Driver
using Constants
using SLutils

export params, kdictGen, genSL, genBZ
# electronic properties
t₁ = 0.1*eV; # In-Bi px-px,py-py
t₂ = 0.4*eV; # In-In 1st NN hopping
t₃ = 0.2*eV; #Bi-Bi 1st NN hopping
t₄ = 0 # obsolete
t₅ = 0.0*eV # In-In 2nd NN hopping
t₆ = 0.0*eV;# Bi-Bi 2nd NN hopping
t₇ = 0.0*eV;# In-Bi further hoppings
t₈ = 0.1*eV;# In-In 2nd nn vertical hoppings
t₉ = 0.0*eV;# In-In px -px, py-py hoppings
t = 1.0*eV; # In-Bi px-px,py-py
ε = 0.0*eV; ε₂ = 0.0*eV #onsite energy for In, Bi

# structural properties
a = 5*Å; c = a
A = [a 0 0; 0 a 0; 0 0 c]
a₁ = A[:,1]; a₂ = A[:,2]; a₃ = A[:,3]
B = transpose(2*π*inv(A))
b₁ = B[:,1]; b₂ = B[:,2]; b₃ = B[:,3]


# modify these in Driver.jl or Runs.jl. Number of supercells in each direction. 
nx = 1; ny = 1; nz = 1; n = nx*ny*nz

function kdictGen(A)
	#B = transpose(π*inv(A))
	# changed for weird def of weyl BZ
	B = transpose(2*π*inv(A))
	#η = (1-a*
	kdict = Dict(
		    "Γ" => B*[ 0  ;    0;   0],
		    "A" => B*[ 1/2;  1/2; 1/2],
		    "M" => B*[ 1/2;  1/2;   0],
		    #"M" => B*[   0;  1/2; 1/2],
		    "Z" => B*[   0;    0; 1/2],
		    "-Z" => B*[  0;    0;-1/2],
		    "X₁" => B*[   1/2;    0; 0],
		    "-X₁" => B*[  -1/2;    0;0],
		    "X₂" => B*[   0;    1/2; 0],
		    "-X₂" => B*[  0;  -1/2;0],
		    "X₃" => B*[   0;    0; 1/2],
		    "-X₃" => B*[  0;    0;-1/2],
		    )
		    #="Γ" => B*[ 0  ;    0;   0],
		    "A" => B*[ 1/2;  1/2; 1/2],
		    "M" => B*[ 1/2;  1/2;   0],
		    "R" => B*[ 0  ;  1/2; 1/2],
		    "X" => B*[ 0  ;  1/2;   0],
		    "-X" => B*[ 0  ;  -1/2;   0],
		    "Z" => B*[ 0  ;    0; 1/2]
		    )=#
		    # or possible for non-tetragonal patterns
		    #"Γ" => B*[ 0  ;    0;   0],
		    #"A" => B*[ 1/2;  1/2; 1/2],
		    #"M" => B*[ 1/2;  1/2;   0],
		    #"R" => B*[ 0  ;  1/2; 1/2],
		    #"X" => B*[ 0  ;  1/2;   0],
		    #"Z" => B*[ 0  ;    0; 1/2],
	    
	return kdict
end

params = (
	  t = t, t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, t₅ = t₅, t₆ = t₆, t₇ = t₇, t₈ = t₈, t₉ = t₉,
	  vf = 10^6, η = 10^-4,
	  ε = ε, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  nx = nx, ny = ny, nz = nz, n = n, norb = 2, nsite = 1,
	  kdict = kdictGen(A)
	  )


function genLatSL(p,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int})
	SLa₁ = SL1[1]*p.a₁ + SL1[2]*p.a₂ + SL1[3]*p.a₃
	SLa₂ = SL2[1]*p.a₁ + SL2[2]*p.a₂ + SL2[3]*p.a₃
	SLa₃ = SL3[1]*p.a₁ + SL3[2]*p.a₂ + SL3[3]*p.a₃
	return SLa₁, SLa₂, SLa₃
end

#=function pruneHoppingType(runtype="")
	if(runtype == "nanopillars")
		return ["z"]
	end
	if(runtype == "domainwall")
		return ["z","y"]
	end
	if(runtype == "device")
		return ["x","y","z"]
	end
	return []
end=#

function genSL(p,nx::Int,ny::Int,nz::Int,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int},runtype::String="",fieldtype="A")
	SLa₁, SLa₂, SLa₃ = genLatSL(p,SL1,SL2,SL3)
	newA = hcat(nx*p.SLa₁,ny*p.SLa₂,nz*p.SLa₃)
	if(nx*ny*nz*p.norb*p.nsite > 64)
		arpack = true
	else
		arpack = false
	end
	SLparams = (
		SLa₁ = newA[:,1], SLa₂ = newA[:,2], SLa₃ = newA[:,3],
		#A = hcat(SLa₁,SLa₂,SLa₃),
		A = newA,
		nx = nx, ny = ny, nz = nz, n=nx*ny*nz,
		kdict = kdictGen(newA),
		runtype=runtype,
		arpack=arpack,
		prune=pruneHoppingType(runtype),
		klist = ["M","Γ","X₁","M","X₂","Γ"],
		fieldtype=fieldtype
	)
	return merge(params,SLparams)
end

function genBZ(p::NamedTuple,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
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
#main(params)

end
