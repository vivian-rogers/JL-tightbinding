push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")

module InBi
using Driver
using Constants

export params, kdictGen
# electronic properties
t₁ = 0.3*eV; # In-Bi px-px,py-py
t₂ = 0.15*eV; # In-In px-py, py-px hopping
t₃ = 0.2*eV; #Bi-Bi px-py, py-px
t₄ = 0 # obsolete
t₅ = 0.03*eV # In-In 2nd nn hopping
t₆ = 0.03*eV;# Bi-Bi 2nd nn hopping
t₇ = 0.02*eV;# In-Bi further hoppings
t₈ = 0.02*eV;# In-In 2nd nn vertical hoppings
t₉ = 0.0*eV;# In-In px -px, py-py hoppings
ε₁ = 0.4*eV; ε₂ = 0.6*eV

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
	  t₁ = t₁, t₂ = t₂, t₃ = t₃, t₄ = t₄, t₅ = t₅, t₆ = t₆, t₇ = t₇, t₈ = t₈, t₉ = t₉,
	  ε₁ = ε₁, ε₂ = ε₂, 
	  a₁ = a₁, a₂ = a₂, a₃ = a₃, A = A, a=a, b=a, c=c,
	  SLa₁ = a₁, SLa₂ = a₂, SLa₃ = a₃,
	  nx = nx, ny = ny, nz = nz, n = n, norb = 2, nsite = 4,
	  kdict = kdictGen(A)
	  )


main(params)

end
