push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./periodicweyl/")
push!(LOAD_PATH, "./src/")

module PWruns

using SparseArrays
using PlotStuff
using LinearAlgebra
using PWHamiltonian
using Bands
using Constants
using UsefulFunctions
using Operators

a = 5*Å

nx = 100; ny = 1; nz = 1;
a₁ = a*[nx;0;0]; a₂ = a*[0;ny;0]; a₃ = a*[0;0;nz]
A = hcat(a₁,a₂,a₃)
B = 2*π*transpose(inv(A))

maxG = Int.(4*round.((1/a)*[maximum(a₁),maximum(a₂),maximum(a₃)])) # highest brillouin zone to sample to. Needs to be big enough to make program not crash.
#sG = 2*maxG+1 # number of G points in PW grid
gcut = 0.5*eV 
Gs, nG = gGrid(B,maxG,gcut) # truncated G grid for hamiltonian

p = (
       arpack = false, η = 10^-8*eV, ng=nG, nG=nG, maxG=maxG, gcut = gcut, A = A, B = B, norb = 2, vf = 10^6, m = 0.1*eV 
         )

sG = (2*maxG[1]+1,2*maxG[2]+1,2*maxG[3]+1)
V = zeros(ComplexF64,sG)
Mx = zeros(ComplexF64,sG)
My = zeros(ComplexF64,sG)
Mz = zeros(ComplexF64,sG)
#=V = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
Mx = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
My = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
Mz = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]);=#

# add stripe domain with periodicity of grid
β = 0.5

# simple helix 
#=
My[gridOffset(p,[1,0,0])] = β*1; My[gridOffset(p,[-1,0,0])] = β*1
Mz[gridOffset(p,[1,0,0])] = β*im; Mz[gridOffset(p,[-1,0,0])] = -β*im
=#


My[gridOffset(p,[1,0,0])] = β*1; My[gridOffset(p,[-1,0,0])] = β*1
My[gridOffset(p,[3,0,0])] = -β*0.1; My[gridOffset(p,[-3,0,0])] = -β*0.1
Mz[gridOffset(p,[1,0,0])] = β*im*0.6; Mz[gridOffset(p,[-1,0,0])] = -β*im*0.6
Mz[gridOffset(p,[3,0,0])] = -β*im*0.2; Mz[gridOffset(p,[-3,0,0])] = β*im*0.2


Mz[gridOffset(p,[0,0,0])] = 0.0000001; 

#Mz[gridOffset(p,[-2,0,0])] = 1/2

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
E, Estates = getBands(klist, kdictGen(A), nk, a, H, p.arpack)
#display(27.2*E)
Q = I(nG)⊗I(p.norb)⊗σ[1]

projStates = expectedValue(Q ,Estates)
plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E)
#plot3DSpectra(H,kdict["M"]+[0;0;0.00*π/a],a,160,160,0.1,0.1)
#plot3DSpectra(H,kdict["M"],a,160,160,0.1,0.1)
println("Done! Press ctrl+d to quit")
end
