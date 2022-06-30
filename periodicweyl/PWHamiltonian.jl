module PWHamiltonian

using LinearAlgebra
using SparseArrays
using UsefulFunctions
using Constants
using Operators

export gtoi, Hβgen, ConstructHamiltonian


function gtoi(p,gvec)
        # gives index in matrix for g vector (-gmax:gmax)
        gmax = p.gcut
        g = gvec[1:3] + gmax
        ng = [gmax[i]*2+1 for i = 1:3]
        return (g[1] + (g[2] + (g[3]*ng[2]))*ng[1]) + 1
end

function itog(p,i)
        
end

function gvectoi(p,gvec)
        # gives index in matrix for g vector (-gmax:gmax), orbital, spin
        gmax = p.gcut
        g = gvec[1:3] + gmax
        ng = [gmax[i]*2+1 for i = 1:3]
        return gvec[5] + (gvec[4] + (g[1] + (g[2] + (g[3]*ng[2]))*ng[1])*p.norb)*2 +1
end

function overlap(p,f::Vector{ComplexF64},gvecs::Vector{Vector{Int}},gvec::Vector{Int},k::Vector{Float64})
        sum = 0
        for g in gvecs
            prod = im
            for ax = 1:3
                θ = 2*π*(g[ax] - gvec[ax])
                #θ = k⋅p.A[:,ax] + 2*π*(g[ax] - gvec[ax])
                #println("θ = $θ")
                if(!(θ≈0))
                        prod *= (exp(-im*θ) - 1)/(θ)
                end
            end
            sum += prod*f[gtoi(p,gvec)]
        end
        return sum
end

function ConstructHamiltonian(p,M)
    Hβ = Hβgen(p,M)
    Hfree = freeElectron(p)
    Hweyl = weylH(p)
    function H(k::Vector{Float64})
        #H =  Hweyl(k)
        #H =  Hweyl(k) .+ I(p.ng*p.norb)⊗σ₂
        H =  Hfree(k) .+ 0.1*I(p.ng*p.norb)⊗σ₃
        #H = Hβ(k) .+ Hfree(k) .+ 2*Diagonal(ones(p.ng*2*p.norb))
        return H
    end
end

function gGrid(p)
        gvecs = Vector{Int}[]
        for gx = -p.gcut[1]:p.gcut[1]
            for gy = -p.gcut[2]:p.gcut[2] 
                for gz = -p.gcut[3]:p.gcut[3]
                    gvec = [gx;gy;gz]
                    #display(gvec)
                    push!(gvecs,gvec)
                end
            end
        end
        return gvecs
end

function weylH(p)
    gvecs = gGrid(p)
    function Hweyl(k::Vector{Float64})
        H = spzeros(ComplexF64,p.ng*p.norb*2,p.ng*p.norb*2)
        # build up a tridiagonal matrix corresponding to the diagonals in the weyl hamiltonian
        # in hilbert space |g>⊗|spin>
        upper = zeros(ComplexF64,p.ng*2-1); diag = zeros(ComplexF64,p.ng*2); lower = zeros(ComplexF64,p.ng*2-1)
        diracterm = Vector{Matrix}(undef,p.ng)
        for ig in eachindex(gvecs)
            G = gvecs[ig]
            gpos = zeros(p.ng); gpos[ig] = 1; gpos = Diagonal(gpos)
            effk = k+p.B*G
            weylAtK = zeros(ComplexF64,2,2)
            for ax = 1:3
                weylAtK .+= effk[ax]*σ[ax]
            end
            H .+= sparse(gpos⊗τ₁⊗weylAtK)
        end
        return (p.vf*(ħ/q)*H.+p.m*τ₃⊗I(p.ng*2))
    end
    return Hweyl
end

function freeElectron(p)
    gvecs = gGrid(p)
    function Hfree(k::Vector{Float64})
        H = zeros(p.ng)
        for g in gvecs
            E = ħ^2*norm(p.B*g + k)^2/(2*m₀*q)
            gi = gtoi(p,g)
            H[gi] = deepcopy(E)
        end
        return Diagonal(H)⊗I(p.norb)⊗I(2)
    end
    return Hfree
end


function Hβgen(p,M)
    gvecs = gGrid(p)
    function Hβ(k::Vector{Float64})
        Hb = spzeros(ComplexF64, p.ng*p.norb*2, p.ng*p.norb*2)
        Mx = M[1]; My = M[2]; Mz = M[3]
        for g in gvecs
            for iorb = 0:(p.norb-1)
                for ispin = 0:1
                    # okay, now take the inner products with M
                    gvec = [g[1];g[2];g[3];iorb;ispin]
                    for mi = 1:3
                        Mᵢ = M[mi]
                        innerprod = overlap(p,Mᵢ,gvecs,g,k)
                        startspin = ispin
                        if(mi < 3) # a little bit of jank
                            startspin = mod(ispin+1,2)
                        end
                        spindict = Dict([0,1]=>1,[0,2]=>-im,[0,3]=>1,[1,1]=>1,[1,2]=>im,[1,3]=>-1)
                        aivec = deepcopy(gvec); aivec[5] = startspin
                        ai = gvectoi(p,aivec)
                        bi = gvectoi(p,gvec)
                        Hb[ai,bi] += innerprod*spindict[[ispin,mi]]
                    end
                end
            end
        end
        return Hb
    end
    return Hβ
end


end
