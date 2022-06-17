push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./periodicweyl/")
push!(LOAD_PATH, "./src/")

module PWruns

using PlotStuff
using LinearAlgebra
using PWHamiltonian
using Bands
using Constants
using UsefulFunctions
using Operators

a = 5*Å

nx = 11; ny = 1; nz = 1;
a₁ = a*[1;0;0]; a₂ = a*[0;1;0]; a₃ = a*[0;0;1]
A = hcat(a₁,a₂,a₃)
B = 2*π*transpose(inv(A))
gcut = [nx-1,ny-1,nz-1]
ng = (gcut[1]*2 + 1)*(gcut[2]*2+1)*(gcut[3]*2+1)

p = (
          gcut=gcut, ng=ng, A = A, B = B, norb = 2, 
         )

Mx = zeros(ComplexF64,ng); Mz = zeros(ComplexF64,ng)
My = zeros(ComplexF64,ng); My[gtoi(p,[1,0,0])] = -im/2; My[gtoi(p,[-1,0,0])] = im/2; My[gtoi(p,[0,0,0])] = 0.2
# ah! a sin wave periodic with the lattice, magnitude 1

function kdictGen(A)
	B = transpose(2*π*inv(A))
	kdict = Dict(
		    "Γ" => B*[ 0  ;    0;   0],
		    "A" => B*[ 1/2;  1/2; 1/2],
		    "M" => B*[ 1/2;  1/2;   0],
		    "Z" => B*[   0;    0; 1/2],
		    "-Z" => B*[  0;    0;-1/2],
		    "X₂" => B*[   0;    1/2; 0],
		    "-X₂" => B*[  0;  -1/2;0],
		    "X₁" => B*[   1/2;    0; 0],
		    "-X₁" => B*[  -1/2;    0;0],
		    )
	return kdict
end

klist = ["M","Γ","X₁","M","X₂","Γ"]

println("Generating periodic field hamiltonian")
H = Hβgen(p,[Mx,My,Mz])

nk = 2^9
println("Getting eigenvalues of 2D weyl lattice between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, kdictGen(A), nk, a, H)
#display(27.2*E)
Q = I(ng)⊗I(p.norb)⊗σ₂

projStates = expectedValue(Q ,Estates)
plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E)
#plot3DSpectra(H,kdict["M"]+[0;0;0.00*π/a],a,160,160,0.1,0.1)
#plot3DSpectra(H,kdict["M"],a,160,160,0.1,0.1)
println("Done! Press ctrl+d to quit")

end
