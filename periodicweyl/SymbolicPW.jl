push!(LOAD_PATH, "../src/")

using SymPy
using Operators
using LinearAlgebra
using BlockDiagonals
using UsefulFunctions

nG = 5; norb = 2
maxG = Int((nG - 1)/2)




a, vf, m, ħ = symbols("a, vf, m, ħ", positive=true)
λ = symbols("λ", integer=true, positive=true)
E = symbols("E")
G = 2*PI/(a*λ)

function weyl()
    kx, ky, kz = symbols("kx, ky, kz")
    k = [kx; ky; kz]
    diag = [τ₁⊗((kx+i*G)*σ₁ .+ ky*σ₂ .+ kz*σ₃) for i = -maxG:maxG]
    #=
    #Rekx, Imkx, ky, kz = symbols("Rekx, Imkx, ky, kz", positive=true)
    k = [Rekx + IM*Imkx; ky; kz]
    diag = [τ₁⊗((Rekx+IM*Imkx+i*G)*σ₁ .+ ky*σ₂ .+ kz*σ₃) for i = -maxG:maxG]
    =#
    #diag = zeros(Sym, nG, nG)
    return Array(BlockDiagonal(diag))
end

H₀ = weyl() .+ m*I(nG)⊗τ₃⊗σ₀ 


display(H₀)
#println("\n\nSpectrum:")
#display(H₀.eigenvals())
println("\n\nCharacteristic polynomial:")
display((E*I(nG*norb*2) .- H₀).det(method="domain-ge"))
#display(H₀.charpoly(E))


