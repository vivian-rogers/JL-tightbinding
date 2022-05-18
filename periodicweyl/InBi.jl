push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module InBi
using Driver
using Constants

export params, kdictGen, genSL
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
ε₁ = 0.0*eV; ε₂ = 0.0*eV #onsite energy for In, Bi

# structural properties
a = 4.8*Å; c = 4.8*Å
A = [a 0 0; 0 a 0; 0 0 c]
a₁ = A[:,1]; a₂ = A[:,2]; a₃ = A[:,3]
B = transpose(2*π*inv(A))
b₁ = B[:,1]; b₂ = B[:,2]; b₃ = B[:,3]


# modify these in Driver.jl or Runs.jl. Number of supercells in each direction. 
nx = 1; ny = 1; nz = 1; n = nx*ny*nz

function kdictGen(A)
	B = transpose(2*π*inv(A))
	kdict = Dict(
		    "Γ" => B*[ 0  ;    0;   0],
		    "A" => B*[ 1/2;  1/2; 1/2],
		    "M" => B*[ 1/2;  1/2;   0],
		    "R" => B*[ 0  ;  1/2; 1/2],
		    "X₁" => B*[ 1/2  ;  0;   0],
		    "X" => B*[ 0  ;  1/2;   0],
		    "Z" => B*[ 0  ;    0; 1/2]
		    )
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
	  t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, t₅ = t₅, t₆ = t₆, t₇ = t₇, t₈ = t₈, t₉ = t₉,
	  vf = 10^7,
	  ε₁ = ε₁, ε₂ = ε₂, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  nx = nx, ny = ny, nz = nz, n = n, norb = 1, nsite = 2,
	  kdict = kdictGen(A)
	  )


function genLatSL(p,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int})
	SLa₁ = SL1[1]*p.a₁ + SL1[2]*p.a₂ + SL1[3]*p.a₃
	SLa₂ = SL2[1]*p.a₁ + SL2[2]*p.a₂ + SL2[3]*p.a₃
	SLa₃ = SL3[1]*p.a₁ + SL3[2]*p.a₂ + SL3[3]*p.a₃
	return SLa₁, SLa₂, SLa₃
end

function genSL(p,nx::Int,ny::Int,nz::Int,SL1::Vector{Int},SL2::Vector{Int},SL3::Vector{Int})
	SLa₁, SLa₂, SLa₃ = genLatSL(p,SL1,SL2,SL3)
	SLparams = (
	SLa₁ = SLa₁, SLa₂ = SLa₂, SLa₃ = SLa₃,
	A = hcat(p.SLa₁,p.SLa₂,p.SLa₃),
	nx = nx, ny = ny, nz = nz, n=nx*ny*nz,
	kdict = kdictGen(A))
	return merge(params,SLparams)
end

#main(params)

end
