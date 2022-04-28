
module PlotStuff
using Plots; pyplot()
using PyPlot
using ColorSchemes
using Constants
using LinearAlgebra


export plotBands, plot2D, SaveFigure, plotVec, plotSurf, plotFunct, plotScatter, plot1D, addLLs, plot3DSpectra

function plotVec(x,yvecs, title)
	nY = size(yvecs)[1]
	for i = 1:nY
		plot(x,yvecs[i])
	end
	title(title)
	gcf()
end

function addLLs(B, fig, nLL=5, Ef=0, vert=true, plot=false)
	
	# yes the /q at the end is weird, use it to convert Joules -> eV
	LLs = [sign(n)*vf*√(2*q*ħ*abs(n)*B)/q for n = -nLL:nLL]
	println("Landau levels: ")
	show(LLs)
	println("\n")
	if(plot)
		figure(fig)
		if(vert)
			for n = 1:(2*nLL+1)
				LL = LLs[n]
				PyPlot.plot([0,10],[LL,LL],linewidth=0.5,c="red")
			end
		else
			for n = 1:(2*nLL+1)
				LL = LLs[n]
				PyPlot.plot([LL,LL],[0,0.05],linewidth=0.5,c="red")
			end
		end
	end
	gcf()
end
		

	

function plot2D(Z,ymin=-4,ymax=4,xlab="",ylab="",name="")
	#nY = size(yvecs)[1]
	#for i = 1:nY
	xmin = 0
	xmax = size(Z)[2]
	
	imshow(Z,cmap="cividis")
	#w, h = PyPlot.figaspect(2.)
	#fig = figure(figsize=(w, h))
    	#ax = fig.add_subplot(111)
    	#ax.set_title(name)
    	#cbar.set_label(zlabel)
    	#ax.set_xlabel(xlab)
    	#ax.set_ylabel(ylab)
    	#fig = PyPlot.imshow(Z, cmap="cividis")
    	#fig = PyPlot.imshow(Z, cmap="cividis",vmin = 0, vmax = 1, origin='lower',interpolation='none')
    	#fig = PyPlot.imshow(Z, cmap="cividis",vmin = 0, vmax = 1, origin="lower",interpolation="none", aspect="auto", extent=[xmin,xmax,ymin,ymax])
    	#divider = PyPlot.make_axes_locatable(ax)
    	#cax = divider.append_axes("right", size="5%", pad=0.05)
    	#PyPlot.colorbar(fig, cax=cax)
	#PyPlot.colorbar(label=name);
	colorbar()
	#=fig, ax = PyPlot.subplots();
	PyPlot.matshow(Arr2D)
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
	#title(title)=#
	gcf()
end

function plotDOS(x,y,xlab="",ylab="")
	#nY = size(yvecs)[1]
	#for i = 1:nY
	fig, ax = PyPlot.subplots();
	PyPlot.plot(x,y)
	#end
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
	#title(title)
	gcf()
end

function plot1D(x,y,xlab="",ylab="")
	#nY = size(yvecs)[1]
	#for i = 1:nY
	fig, ax = PyPlot.subplots();
	PyPlot.plot(x,y)
	#end
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
	#title(title)
	gcf()
end



function plot3DSpectra(H,k₀,a,nx,ny,dx,dy)
	nE = size(H([0,0,0]))[1]
	# generate k grid
	k₀x = k₀[1]; k₀y = k₀[2]; k₀z = k₀[3]
	x = range(k₀x - (π/(a))*dx/2,stop= k₀x + (π/(a))*dx/2,length=nx)
	y = range(k₀y - (π/(a))*dy/2,stop= k₀y + (π/(a))*dy/2,length=ny)
	#y = range(k₀y - dy/2,stop= k₀y + dy/2,length=nx)

	println("Plotting 3D spectra around k₀")
	#fig = Plots.scatter(x,y,f, camera=(-30,30))
	function f(kx,ky)
		E = eigvals(H([kx,ky,k₀z]))
		return minimum(abs.([E[1]-E[3], E[2]-E[3], E[2]-E[4]]))
		#return minimum(abs.([E[1]-E[2], E[2]-E[3], E[3]-E[4]]))
	end
	#plot
	#Plots.plot(x,y,f,st=:surface, camera=(10,10), palette = :balance, alpha = 0.3)
	#Plots.plot(x,y,f,st=:surface, camera=(10,10), palette = :balance, alpha = 0.3, xlim = (minimum(x),maximum(x)), ylim = (minimum(y),maximum(y)))
	Plots.heatmap!(x,y,f, palette = :balance, xlim = (minimum(x),maximum(x)), ylim = (minimum(y),maximum(y)), xlab="k₁ (π/a)", ylab="k₂ (π/a)",title="ΔE between band 2 and 3", clim =(0,0.02))
	#Plots.heatmap!(x,y,f, palette = :balance, alpha = 0.3, xlim = (minimum(x),maximum(x)), ylim = (minimum(y),maximum(y)))
	Plots.plot!([0,π/a,π/a],[π/a,π/a,0], lw=1)
	#Plots.plot([0,π/a,π/a],[π/a,π/a,0],[0,0,0])
	#=for iE = 1:nE
		f(kx,ky) = (eigvals(H([kx,ky,k₀z])))[iE]
		println("eigvl #$iE...")
		Plots.plot!(x,y,f,st=:surface, camera=(10,10), palette = :balance, alpha = 0.3)
	end=#
	#gcf()
	gui()
	#display(fig)
	#display(fig)
	#spectra = map(k->eigvals(H(k)), kgrid)
	#for iE = 1:nE
		
		#for ik in spectra:
			
end


function plotBands(klist, nk, E, projStates, name="")
	#Plots.pyplot()
	nSymPts = size(klist)[1]
	indices = LinRange(0,nSymPts-1, size(E)[1])
	nE = size(E)[2]
	fig, ax = PyPlot.subplots();
	#display(E[:,1])
	#display(plot!(indices,E[1,:]))
	#display(plot!(indices,E[:,2]))
	#Eplot = transpose(E)
	kSymPts = [i for i =0:(nSymPts-1)]
	for kTick in kSymPts
		PyPlot.plot([kTick,kTick],[-30,30],c="#666666",lw=0.5)
	end
	Pyplot.xticks(kSymPts,klist)
	xlim(0, nSymPts-1)
	maxE = maximum(E)
	minE = minimum(E)
	ylabel("E - Ef (eV)")
	ylim(minE+0.1,maxE-0.1)


	#set_cmap("rainbow")
	for iE = 1:nE
		Evals = collect(E[:,iE])
		Projvals = collect(projStates[:,iE])
		#Projvals = collect(projStates[:,iE])
		#scatter(indices,Evals,c=Projvals, s=0.9)
		PyPlot.scatter(indices,Evals,c=Projvals, cmap="coolwarm", s=0.9, vmin = -1, vmax = 1)
		#scatter(indices,Evals,c=Projvals, vmin=0, vmax=1,s=0.9)
		#display(plot!(indices,Evals))
	end
	PyPlot.colorbar(label=name);
	gcf()
end

function plotPoints(Rvals,z,xlab="",ylab="",zlab="",name="",cmap= :redgreensplit,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	surf = scatter(Rvals[:,1],Rvals[:,2],z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	gui(surf)
end

function plotSurf(x,y,z,xlab="",ylab="",name="",cmap= :redgreensplit,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	#heatmap!
	gui(surf)
end

#function nameGen

function SaveFigure(fig,path,name="",type=".png")
	fig.savefig(path*"/"*name*type)
	close(fig)
end

function plotScatter(R,z,xlab="",ylab="",name="",cmap="inferno",xyscale=(1/nm),zscale=1,save=false)
	dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	C = 500
	width = C
	height = C*dy/dx
	fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=Tuple([zscale*zi for zi in z]));
	PyPlot.xlabel(xlab);
	#fig.colorbar(zplot, ax=ax);
	PyPlot.colorbar(zplot, label=name);
	PyPlot.ylabel(ylab);
	gcf()
	#show(surf);
end


function plotFunct(R,f::Function,xlab="",ylab="",name="",cmap="bwr",xyscale=(1/nm),zscale=1,save=false)
	dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	C = 500
	width = C
	height = C*dy/dx
	z = f.(R)
	#show(z)
	#surf = Plots.scatter(R[:][1],R[:][2],z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	
	#surf = Plots.scatter(R[:,1],R[:,2], marker=:c, color=cmap, zcolor=z, size=(width,height))
	#show([r[1] for r in R])
	#surf = Plots.scatter([r[1] for r in R],[r[2] for r in R], marker=:c, color=cmap, zcolor=z, size=(width,height))
	fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), cmap=cmap, c=zscale*z);
	PyPlot.xlabel(xlab);
	#fig.colorbar(zplot, ax=ax);
	PyPlot.colorbar(zplot, label=name);
	PyPlot.ylabel(ylab);
	gcf()
end


function plotBands(klist, nk, E)
	#Plots.pyplot()
	nSymPts = size(klist)[1]
	indices = LinRange(0,nSymPts-1, size(E)[1])
	nE = size(E)[2]
	#display(E[:,1])
	#display(plot!(indices,E[1,:]))
	#display(plot!(indices,E[:,2]))
	#Eplot = transpose(E)
	kSymPts = [i for i =0:(nSymPts-1)]
	for kTick in kSymPts
		PyPlot.plot([kTick,kTick],[-30,30],c="#666666",lw=0.5)
	end
	PyPlot.xticks(kSymPts,klist)
	xlim(0, nSymPts-1)
	maxE = maximum(E)
	minE = minimum(E)
	#maxE = Emax
	#minE = -Emax
	ylabel("E - Ef (eV)")
	ylim(minE,maxE)
	for iE = 1:nE
		Evals = collect(E[:,iE])
		#plot(indices,Evals)
		PyPlot.scatter(indices,Evals,s=1)
		#display(plot!(indices,Evals))
	end
	gcf()
end

end