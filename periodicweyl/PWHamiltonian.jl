module PWHamiltonian

using LinearAlgebra
using SparseArrays

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
                θ = k⋅p.A[:,ax] + 2*π*(g[ax] - gvec[ax])
                #println("θ = $θ")
                if(!(θ≈0))
                        prod *= im*(exp(-im*θ) - 1)/(θ)
                end
            end
            sum += prod*f[gtoi(p,gvec)]
        end
        return sum
end

function ConstructHamiltonian(p,M)
    Hβ = Hβgen(p,M)
    function H(k::Vector{Float64})
        Hᵦ = Hβ(k)
        return Hᵦ
    end
end

function gGrid(p)
        gvecs = Vector{Int}[]
        for gx = -p.gcut[1]:p.gcut[1]
            for gy = -p.gcut[2]:p.gcut[2] 
                for gz = -p.gcut[3]:p.gcut[3]
                    gvec = [gx;gy;gz]
                    display(gvec)
                    push!(gvecs,gvec)
                end
            end
        end
        return gvecs
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
                        aivec = deepcopy(gvec); aivec[5] = startspin
                        ai = gvectoi(p,aivec)
                        bi = gvectoi(p,gvec)
                        Hb[ai,bi] += innerprod
                    end
                end
            end
        end
        return Hb
    end
    return Hβ
end


end
