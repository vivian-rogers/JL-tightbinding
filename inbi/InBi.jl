push!(LOAD_PATH, "./")
push!(LOAD_PATH, "./src/")
push!(LOAD_PATH, "../src/")

module InBi
using Driver
using Constants

export params, kdictGen
# electronic properties
t₁ = 2.5*eV; t₂ = 1.0*eV; t₃ = 1.5*eV; t₄ = 0.6*eV
ε₁ = -0.2*eV; ε₂ = 0.4*eV

# structural properties
a = 4.8*Å; c = 4.8*Å
A = [a 0 0; 0 a 0; 0 0 c]
a₁ = A[:,1]; a₂ = A[:,2]; a₃ = A[:,3]
B = transpose(2*π*inv(A))
b₁ = B[:,1]; b₂ = B[:,2]; b₃ = B[:,3]

nx = 1; ny = 1; nz = 1; n = nx*ny*nz

function kdictGen(A)
	B = transpose(2*π*inv(A))
	kdict = Dict(
		    "Γ" => B*[ 0  ;    0;   0],
		    "A" => B*[ 1/2;  1/2; 1/2],
		    "M" => B*[ 1/2;  1/2;   0],
		    "R" => B*[ 0  ;  1/2; 1/2],
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
	  t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, ε₁ = ε₁, ε₂ = ε₂, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  nx = nx, ny = ny, nz = nz, n = n, norb = 2, nsite = 4,
	  kdict = kdictGen(A)
	  )


main(params)

end
