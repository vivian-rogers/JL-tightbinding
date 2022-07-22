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

nx = 1; ny = 1; nz = 1;
a₁ = a*[nx;0;0]; a₂ = a*[0;ny;0]; a₃ = a*[0;0;nz]
A = hcat(a₁,a₂,a₃)
B = 2*π*transpose(inv(A))

maxG = 50 # highest brillouin zone to sample to
sG = 2*maxG+1 # number of G points in PW grid
gcut = 5.0*eV
Gs, nG = gGrid(B,maxG,gcut) # truncated G grid for hamiltonian

p = (
       ng=nG, nG=nG, maxG=maxG, gcut = gcut, A = A, B = B, norb = 1, vf = 10^6, m = 0.5 
         )

V = zeros(ComplexF64,sG,sG,sG); Mx = zeros(ComplexF64,sG,sG,sG); Mz = zeros(ComplexF64,sG,sG,sG); My = zeros(ComplexF64,sG,sG,sG); 

# add cosine with periodicity of grid
#Mz[gridOffset(p,[0,0,0])] = 0.5;
Mz[gridOffset(p,[0,0,0])] = 0.5;
#Mx[gridOffset(p,[0,0,0])] = 0.5;
#β = 0.5
#V[gridOffset(p,[1,0,0])] = β*1; V[gridOffset(p,[-1,0,0])] = β*1
#My[gridOffset(p,[1,0,0])] = β*1; My[gridOffset(p,[-1,0,0])] = β*1
#Mz[gridOffset(p,[1,0,0])] = β*im; Mz[gridOffset(p,[-1,0,0])] = -β*im
#Mz[gridOffset(p,[2,0,0])] = 1/2; Mz[gridOffset(p,[-2,0,0])] = 1/2

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
		    "X₃" => B*[  0;    0;1/2],
		    )
	return kdict
end

klist = ["M","Γ","X₁","M","X₂","Γ","X₃"]

println("Generating periodic field hamiltonian")
H = ConstructHamiltonian(p,[V,Mx,My,Mz])
#H = Hβgen(p,[Mx,My,Mz])

nk = 2^8
println("Getting eigenvalues of 2D weyl lattice between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, kdictGen(A), nk, a, H)
#display(27.2*E)
Q = I(nG)⊗I(p.norb)⊗σ[3]

projStates = expectedValue(Q ,Estates)
plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E)
#plot3DSpectra(H,kdict["M"]+[0;0;0.00*π/a],a,160,160,0.1,0.1)
#plot3DSpectra(H,kdict["M"],a,160,160,0.1,0.1)
println("Done! Press ctrl+d to quit")

end
