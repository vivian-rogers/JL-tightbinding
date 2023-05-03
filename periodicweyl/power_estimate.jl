
nm = 10^-9; cm = 10^-2; pico = 10^-12; giga = 10^9; nano = 10^-9; mega = 10^6

function performance(p)
	#MR = MR(p)
	Δt1 = Δt_1state(p)
	Δtr = Δt_reset(p)
	nstates = p.L/p.d + 1
	println("Performance of " * p.name)
	println("Magnetoresistance = $(MR(p)) %, ΔE_reset = $(ΔE(p)/pico) pJ")
	println("Time for full reset = $(Δtr/nano) ns and reset freq = $(Δtr^-1/mega) MHz")
	println("Time for moving 1 state = $(Δt1/nano) ns and reset freq = $(Δt1^-1/mega) MHz, nstates = $(nstates)")
end

function MR(p)
	return (p.L/p.d+1)*(p.Rdw/p.R₀)*100
end

function Δt_1state(p)
	return p.d/p.v
end

function Δt_reset(p)
	return p.L/p.v
end

function ΔE(p)
	Rhi = (p.L/p.d + 1)*p.Rdw+p.R₀
	return (p.L/p.v)*(p.J*p.t*p.W)^2*Rhi
end

function genNewParams(p, name, t, W, L, d, Rdw)
	ρ₀ = p.W*p.t*p.R₀/p.L
	ρDW = p.W*p.t*Rdw
	return merge(p, (L = L, t =t, W= W, d=d, name = name, R₀ = L*ρ₀/(t*W), Rdw = ρDW/(W*t)))
end
function genNewParams(p, name, t, W, L)
	ρ₀ = p.W*p.t*p.R₀/p.L
	ρDW = p.W*p.t*p.Rdw
	return merge(p, (L = L, t =t, W= W, name = name, R₀ = L*ρ₀/(t*W), Rdw = ρDW/(W*t)))
end

#pMn3Sn_best = (name="Mn3Sn best case", L = 35*1000*nm, v = 1.5, d=500*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 1000*nm, Rdw = 6, R₀ =309)
pMn3Sn_reported = (name="Mn3Sn reported params", L = 35*1000*nm, v = 1.5, d=800*nm, J = 1.5*10^5*cm^-2, t = 500*nm, W = 1000*nm, Rdw = 6, R₀ =309)
pMn3Sn_real = genNewParams(pMn3Sn_reported,"Mn3Sn_real", 40*nm, 1000*nm, 10000*nm)
pMn3Sn_best = genNewParams(pMn3Sn_reported,"Mn3Sn_best", 40*nm, 300*nm, 5000*nm, 300*nm, 15)
#p = (name="Mn3Sn best case", L = 35*1000*nm, v = 1.5, d=500*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 1000*nm, Rdw = 6, R₀ =309)
#p = (L = 35*1000*nm, v = 1.5, d=800*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 1000*nm, Rdw = 6, R₀ =309)
pCo3Sn2S2_reported = (name="Co3Sn2S2_reported", L = 600*nm, v = 1.5, d=80*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 200*nm, Rdw = 2, R₀ =53)
pCo3Sn2S2_real = genNewParams(pCo3Sn2S2_reported, "Co3Sn2S2_realistic",40*nm,500*nm, 2*1000*nm)
pCo3Sn2S2_best = genNewParams(pCo3Sn2S2_reported, "Co3Sn2S2 best case",40*nm,500*nm, 5*1000*nm, 40*nm, 5)
#pCo3Sn2S2_real = (name="Co3Sn2S2_real",L = 3000*nm, v = 1.5, d=80*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 500*nm, Rdw = 2, R₀ =53)
#pCo3Sn2S2_best = (name="Co3Sn2S2_best",L = 3000*nm, v = 4.5, d=80*nm, J = 1.5*10^5*cm^-2, t = 40*nm, W = 500*nm, Rdw = 5, R₀ =53)
#p = (L = 600*nm, v = 1.5, d=2*nm, J = 1.5*10^5*cm^-2, t = 10*nm, W = 100*nm, Rdw = 0.5, R₀ =50)
for p ∈ [pMn3Sn_reported, pMn3Sn_real, pMn3Sn_best, pCo3Sn2S2_reported, pCo3Sn2S2_real, pCo3Sn2S2_best]
	performance(p)
end
