push!(LOAD_PATH, "./src/")

#using Plots
using nearestNeighbors
using Nanospheres
using Constants
using UsefulFunctions
using PlotStuff
using GMaterialParameters

λ = 84; σ = 3*10^-9; Δh0 = 2*10^-9; fᵤ = 1

#x,y,z,fz = NSprofileGen(λ,Δh0,σ)
#R = MtoA(Rvals(λ))
#dux_dx = testhess(fz)
n = 5
klist, kweights = kgrid(n,λ)
plotScatter(klist,kweights)
#psB, u = strain(λ,Δh0,fᵤ,σ,fz)
#pseudoB, dux_dx, ε_eff, dux_dx, dh_dx = strain(λ,Δh0,fᵤ,σ,fz)
#testh = fz.(R)
#show(R)
#testf = dux_dx.(R)
#plotFunct(R,pseudoB,"x position (nm)","y position (nm)","Bₚₛ (T)")
#plotFunct(R,fz,"x position (nm)","y position (nm)","height (nm)","inferno",(1/nm),(1/nm))
#plotFunct(R,ε_eff,"x position (nm)","y position (nm)","eff. ε")
#plotFunct(R,dux_dx,"x position (nm)","y position (nm)","dux_dx")
#plotFunct(R,dh_dx,"x position (nm)","y position (nm)","dh_dx")
#plotFunct(R,pseudoB,"x position (nm)","y position (nm)","Bₚₛ (T)")
#plotFunct(R,pseudoB,"x position (nm)","y position (nm)","pseudoB","bwr",(1/nm),1)
#plotFunct(R,pseudoB,"x position (nm)","y position (nm)","pseudoB",:inferno,(1/nm),(1/nm))

#pseudoB = 
#surf = surface(x,y,z)

#plotSurf((1/nm)*x,(1/nm)*y,(1/nm)*z,"x position (nm)", "y position (nm)","Height profile", :plasma)
#surf = plot(x,y,z,st=:surface)
#gui(surf)
#surface(x,y,z)

