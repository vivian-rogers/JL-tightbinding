push!(LOAD_PATH, "./")
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./periodicweyl/")
push!(LOAD_PATH, "./src/")

#module PWruns

using SparseArrays
using PlotStuff
using LinearAlgebra
using PWHamiltonian
using Bands
using Constants
using UsefulFunctions
using Operators
using SLutils
using DOS

a = 5*Å

nx = 100; ny = 1; nz = 1;


function paramGen(nx::Int,ny::Int,nz::Int)
    a₁ = a*[nx;0;0]; a₂ = a*[0;ny;0]; a₃ = a*[0;0;nz]
    A = hcat(a₁,a₂,a₃)
    B = 2*π*transpose(inv(A))

    maxG = Int.(round.((4/a)*[maximum(a₁),maximum(a₂),maximum(a₃)])) # highest brillouin zone to sample to. Needs to be big enough to make program not crash.
    #sG = 2*maxG+1 # number of G points in PW grid
    #gcut = 8.5*eV 
    gcut = [5;0;0]
    Gs, nG = gGrid(B,maxG,gcut) # truncated G grid for hamiltonian

    sG = (2*maxG[1]+1,2*maxG[2]+1,2*maxG[3]+1)
    params = (
           arpack = false, η = 10^-3*eV, μ_disorder=0.001*eV, ng=nG, nG=nG, sG = sG,
           maxG=maxG, gcut = gcut, A = A, B = B, norb = 2, vf = 10^6, m = 0.5*eV, t = 1*eV
             )
    return params
end

p = paramGen(nx,ny,nz)

V = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
Mx = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
My = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
Mz = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
#Mx = spzeros(ComplexF64,p.sG)
#My = spzeros(ComplexF64,p.sG)
#Mz = spzeros(ComplexF64,p.sG)
#=V = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
Mx = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
My = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]); 
Mz = zeros(ComplexF64,maxG[1],maxG[2],maxG[3]);=#

# add stripe domain with periodicity of grid
β = 1.0
#V[gridOffset(p,[1,0,0])] = β*1; V[gridOffset(p,[-1,0,0])] = β*1


#simple neel helix 
#Mx[gridOffset(p,[1,0,0])] = β*1; Mx[gridOffset(p,[-1,0,0])] = β*1
#My[gridOffset(p,[1,0,0])] = β*im; My[gridOffset(p,[-1,0,0])] = -β*im

# simple bloch helix 

My[gridOffset(p,[1,0,0])] = β/2; My[gridOffset(p,[-1,0,0])] = β*1/2
Mz[gridOffset(p,[1,0,0])] = β*im/2; Mz[gridOffset(p,[-1,0,0])] = -β*im/2


#simplest bloch lattice

#=
My[gridOffset(p,[1,0,0])] = β*1; My[gridOffset(p,[-1,0,0])] = β*1
My[gridOffset(p,[3,0,0])] = -β*0.1; My[gridOffset(p,[-3,0,0])] = -β*0.1
Mz[gridOffset(p,[1,0,0])] = β*im*0.6; Mz[gridOffset(p,[-1,0,0])] = -β*im*0.6
Mz[gridOffset(p,[3,0,0])] = -β*im*0.2; Mz[gridOffset(p,[-3,0,0])] = β*im*0.2
=#


#Mx[gridOffset(p,[1,0,0])] = 1/2; Mx[gridOffset(p,[-1,0,0])] = 1/2


Mx[gridOffset(p,[0,0,0])] = 0.001; 
#Mz[gridOffset(p,[0,0,0])] = 1.0; 

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

klist = ["M","-X₁","Γ","X₁","M","X₂","Γ","X₃"]

println("Generating periodic field hamiltonian")
H = ConstructHamiltonian(p,[V,Mx,My,Mz])
#H = Hβgen(p,[Mx,My,Mz])

nk = 2^8
println("Getting eigenvalues of 2D weyl lattice between k = ")
show(klist)
println("...")

#=
Q = I(p.nG)⊗σ[1]⊗I(2)
E, Estates = getBands(klist, kdictGen(p.A), nk, a, H, p.arpack)

projStates = expectedValue(Q ,Estates)
plotBands(klist,nk,E, projStates)
=#

γ⁵ = σ[1]⊗I(2)
γᴸ= (1/2)*(I(4) .- γ⁵)
nk = 25
#kslice(p,H,0.24,"x",nk,nk,1.0)
#energySurface(p,H,0.20,3,nk,nk)
neigs = 8*p.arpack + p.nG*p.norb*2*(!p.arpack)

eigSurface(p,H,I(p.nG*p.norb)⊗σ[3],neigs,"y",nk,nk,0.0)
#eigSurface(p,ConstructHamiltonian(p,[V,0.55*Mx,0.25*My,0.25*Mz]),I(p.nG)⊗γ⁵,neigs,"z",nk,nk,0.0)

#complexEnergySurface(p,H,0.0,400,nk,nk,2,3)
#complexEigSurface(p,ConstructHamiltonian(p,[V,Mx,My,Mz]),I(p.nG)⊗I(p.norb)⊗σ[2],neigs,2,nk,nk,2*nm,[0.0;0.0;0.0])
#eigSurface(p,H,I(p.nG)⊗γ⁵,neigs,"z",4,nk,0.0)
#eigSurface(p,H,I(p.nG)⊗I(p.norb)⊗σ[2],neigs,"z",7,nk,0.0)


function arbβToH(β::Float64, λ::Int)
    V = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    Mx = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    My = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    Mz = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    β = 1.0
    My[gridOffset(p,[1,0,0])] = β/2; My[gridOffset(p,[-1,0,0])] = β*1/2
    Mz[gridOffset(p,[1,0,0])] = β*im/2; Mz[gridOffset(p,[-1,0,0])] = -β*im/2
    p = paramGen(λ,1,1)
    H = ConstructHamiltonian(p,[V,Mx,β*My,β*Mz])
    return H
end

function arbβToH(β::Float64, λ::Int)
    V = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    Mx = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    My = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    Mz = zeros(ComplexF64,p.sG[1],p.sG[2],p.sG[3])
    β = 1.0
    My[gridOffset(p,[1,0,0])] = β/2; My[gridOffset(p,[-1,0,0])] = β*1/2
    Mz[gridOffset(p,[1,0,0])] = β*im/2; Mz[gridOffset(p,[-1,0,0])] = -β*im/2
    p = paramGen(λ,1,1)
    H = ConstructHamiltonian(p,[V,Mx,β*My,β*Mz])
    return H
end

function βconstβToH(βconst::Float64, βhelix::Float64)
#function λβToH(λ::Int, β::Float64)
    #p = paramGen(λ,1,1)
    H = ConstructHamiltonian(p,[V,βconst*Mx,βhelix*My,βhelix*Mz])
    return H
end

#Sweep2DSurf(findBandgap(p),βconstβToH, [βconst for βconst = 0.0:0.05:1.0], [βhelix for βhelix = 0.0:0.05:1.0],   "β x constant (eV)", "β y,z helix (eV)", "Bandgap (eV)")

#Sweep2DSurf(findBandgap(p),βconstβToH, [βconst for βconst = 0.0:0.05:1.0], [βcos for βcos = 0.0:0.05:1.0],   "β x constant (eV)", "β y,z helix (eV)", "Bandgap (eV)")
#Sweep2DSurf(bandgap,λβToH, [βconst for βconst = 0.0:0.1:1.0], [λ for λ = 5:25:100], "√(M₂² + M₃²) (eV)", "λ")
#energySurface(p,H,0.25,1,nk,nk)
#plotBands(klist,nk,E, projStates)
#plotBands(klist,nk,E)
println("Done! Press ctrl+d to quit")
