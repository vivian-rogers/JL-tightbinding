push!(LOAD_PATH, "../src/")

using AbstractAlgebra
using Symbolics
using Operators
using LinearAlgebra
using BlockDiagonals
using UsefulFunctions
using Constants

nG = 3; norb = 2
maxG = Int((nG - 1)/2)
nλ = 100; na = 5*Å;
fG = 2*π/(na*nλ)




(a, vf, m, ħ, β, λ, E) = @variables a::Real vf::Real m::Real ħ::Real β::Real λ::Real E::Complex
G = 2*π/(a*λ)

# generate the weyl hamiltonian
kx, ky, kz = @variables kx::Complex ky::Real kz::Real
diag = [τ₁⊗((kx+i*G)*σ₁ .+ ky*σ₂ .+ kz*σ₃) for i = -maxG:maxG]
Hweyl =  Array(BlockDiagonal(diag))


# define the periodic exchange field
function gtoi(g::Int)
    return Int(maxG+1+g)
end
Mx = zeros(Complex,nG)
My = zeros(Complex,nG)
Mz = zeros(Complex,nG)
#My[gtoi(-1)] = im*β/2;    My[gtoi(1)] = -im*β/2
Mz[gtoi(-1)] = β/2; Mz[gtoi(1)] = β/2
M = [Mx,My,Mz]

# generate the periodic exchange field term
Hᵦ = zeros(Complex,nG*norb*2,nG*norb*2)
for ax = 1:3
    Mi = M[ax]
    Mᵢ(index) = (if(index < 1 || index > (2*maxG+1)) return 0 else return Mi[index] end)
    Gterm = zeros(Sym,nG,nG)
    for gb = -maxG:maxG
        for ga = -maxG:maxG
            ia = gtoi(ga); ib = gtoi(gb);
            Gterm[ib,ia] = Mᵢ(gtoi(gb-ga))
        end
    end
    Hᵦᵢ = Gterm⊗Num.(τ₀⊗σ[ax])
    #display(Hᵦᵢ)
    Hᵦ .+= Hᵦᵢ
end

display(Hᵦ)
println("\nHᵦ ^\n")

H₀ = vf*Hweyl .+ m*I(nG)⊗τ₃⊗σ₀ .+ Hᵦ



#Rekx, Imkx, ky, kz = symbols("Rekx, Imkx, ky, kz", positive=true)
#diag = [τ₁⊗((Rekx+IM*Imkx+i*G)*σ₁ .+ ky*σ₂ .+ kz*σ₃) for i = -maxG:maxG]



display(H₀)
#println("\n\nSpectrum:")
#display(H₀.eigenvals())
println("\n\nGenerating characteristic polynomial...")
#@time charpoly = (E*I(nG*norb*2) .- H₀).det(method="berkowitz")
#=
@time charpoly = (E*I(nG*norb*2) .- H₀).det(method="LU")
#atΓ = charpoly.subs(Rekx,0)

#display(charpoly)
atΓ = charpoly.subs([(kz,0)])
#atΓ = charpoly.subs([(kz,0),(E,0)])
# hartree unit definitions
#atΓ = atΓ.subs([(ħ,hbar),(vf,10^6*MpS),(m,0.5*eV),(λ,nλ),(a,10*a₀),(β,1*eV)])
atΓ = simplify(atΓ)
atΓ = cancel(atΓ)
atΓ = collect(atΓ,kx)
display(atΓ)
atΓ = atΓ.subs(ky,0)
polyΓ = (atΓ.evalf()).as_poly(domain="C")
coeffs = SymPy.all_coeffs(polyΓ)

#atΓ = cancel(atΓ)
#atΓ = simplify(charpoly.subs([(ky,0),(kz,0),(E,0.05)]))
#display(atΓ)
#atΓ = atΓ.subs([(ħ,1),(vf,1),(m,0.5),(λ,nλ),(a,1),(β,0.00000000000001)])
#decay = atΓ.subs(ħ,Constants.ħ/Constants.q)
#decay = decay.subs(vf,10^6)
#decay = decay.subs(m,0.5)
println("\n\nCharacteristic polynomial @ Γ:")
#println("Too hard to display!")
#atΓ = atΓ.evalf()
#=
display(atΓ)
println("\n\nSolving...")
#sols = N.(solve(Eq(atΓ,0), warn=true, rational=false, implicit=true, particular=true))
#sols = nsolve(atΓ,IM)

polyΓ = ((atΓ).evalf()).as_poly(domain="C")
hbar = 1; q = 1; a₀ = 1; nm = 18.897*a₀; Eh = 1; eV =Eh/27.211; MpS = (a₀*Eh/hbar)/(2.18769*10^6)
#atΓ = atΓ.subs([(vf, 10^6*MpS),(m, 0.5*eV),(λ,nλ),(a,0.5*nm),(β,1*eV),(ky,0),(E,0)])
polyΓ = polyΓ.subs([(vf, 10^6*MpS),(m, 0.5*eV),(λ,nλ),(a,0.5*nm),(β,1*eV),(ky,0),(E,0)])
sols = N.(solve(Eq(polyΓ,0), warn=true, rational=false, implicit=true, particular=true))
display(sols)
sols = roots(polyΓ)
sols = [key for (key, value) in sols]
#sols = roots(polyΓ,multiple=true)
#sols = nsolve(atΓ,im*(2π/na)+nG*0.1)
#sols = linsolve((atΓ,0),kx)
#sols = (linsolve(Eq(atΓ,0),rational=false))
#sols = ComplexF64.(sols)
display(sols)
#realk = mod.(real.(sols)./(fG/2),2)
=#
realk = real.(sols)./((fG)/2)
imagk = inv.(nm*imag.(sols))
println("\n\nRe(sols) = $realk * 𝐗 ")
println("\nIm(sols) -> decay rates = $imagk nm")


