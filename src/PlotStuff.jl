
module PlotStuff
using Plots; pyplot()
using PyPlot
using ColorSchemes
using Constants
using LinearAlgebra
using Distributed
using ProgressBars
using Measures

export plotBands, plot2D, SaveFigure, SavePlots, plotVec, plotSurf, plotFunct, plotScatter, plot1D, addLLs, plot3DSpectra, plot3Dvectors, plotMat


fsize=15; fname="arial"
rc("font", family=fname, size=fsize)
pyplot()
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
		

	
function plotMat(Z,xlab="",ylab="",name="")
	#nY = size(yvecs)[1]
	#for i = 1:nY
        pyplot()
        h = heatmap(Z,cmap="inferno",xlabel=xlab,ylabel=ylab)
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
	#colorbar()
	#=fig, ax = PyPlot.subplots();
	PyPlot.matshow(Arr2D)
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
	#title(title)=#
	gui(h)
end

function plot2D(Z,ymin=-4,ymax=4,xlab="",ylab="",name="")
	#nY = size(yvecs)[1]
	#for i = 1:nY
        pyplot()
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
	#colorbar()
	#=fig, ax = PyPlot.subplots();
	PyPlot.matshow(Arr2D)
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
	#title(title)=#
	gui()
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

function plot1D(x::Vector,ys::Vector{Vector{Float64}},xlab::String="",ylab::String="",legend::Vector=[],yaxis=:log)
	#nY = size(yvecs)[1]
	#for i = 1:nY
        pyplot()
        for i in eachindex(ys)
            y = ys[i]
            if(i==1)
                Plots.plot(x,y,xlabel=xlab,ylabel=ylab,label=legend[i],yscale=yaxis,color=i, margin = 15mm)
                Plots.scatter!(x,y,xlabel=xlab,ylabel=ylab,label=false,yscale=yaxis,color=i)
            else
                Plots.plot!(x,y,xlabel=xlab,ylabel=ylab,label=legend[i],color=i,yscale=yaxis)
                Plots.scatter!(x,y,xlabel=xlab,ylabel=ylab,label=false,color=i,yscale=yaxis)
            end
        end
        #    #if(legend == [])
        #    Plots.plot(x,ys,xlabel=xlab,ylabel=ylab,label=legend)
        #else
        #    Plots.plot(x,ys,xlabel=xlab,ylabel=ylab)
        #end
        xmin = minimum(minimum.(x));
       xmax = maximum(maximum.(x));
       ymin = minimum(minimum.(ys));
       ymax = maximum(maximum.(ys));
        ylim(ymin,ymax)
        xlim(xmin,xmax)
        gui()
        #fig, ax = PyPlot.subplots();
	#PyPlot.plot(x,y)
	#end
	#PyPlot.ylabel(ylab);
	#PyPlot.xlabel(xlab);
        #PyPlot.xlim(xmin,xmax)
        #PyPlot.ylim(ymin,ymax)
	#gcf()
end

function plot1D(x,y,xlab="",ylab="",slab=false,xmin::Float64=0.0,xmax::Float64=1.0,ymin::Float64=-1.0,ymax::Float64=1.0)
	#nY = size(yvecs)[1]
	#for i = 1:nY
        f = Plots.plot(x,y,xlabel=xlab,ylabel=ylab, guidefont=font(fsize,fname),label=slab)
        ylim(ymin,ymax)
        xlim(xmin,xmax)
        gui()
        return f
        #fig, ax = PyPlot.subplots();
	#PyPlot.plot(x,y)
	#end
	#PyPlot.ylabel(ylab);
	#PyPlot.xlabel(xlab);
        #PyPlot.xlim(xmin,xmax)
        #PyPlot.ylim(ymin,ymax)
	#gcf()
end

function pyplot1D(x,y,xlab="",ylab="",xmin::Float64=0.0,xmax::Float64=1.0,ymin::Float64=-1.0,ymax::Float64=1.0)
	#nY = size(yvecs)[1]
	#for i = 1:nY
	fig, ax = PyPlot.subplots();
	PyPlot.plot(x,y)
	#end
	PyPlot.ylabel(ylab);
	PyPlot.xlabel(xlab);
        PyPlot.xlim(xmin,xmax)
        PyPlot.ylim(ymin,ymax)
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
	PyPlot.xticks(kSymPts,klist)
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
		#PyPlot.scatter(indices,Evals,c=Projvals, cmap="coolwarm", s=0.9)
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
	height = C
        height = C*dy/dx
	surf = scatter(Rvals[:,1],Rvals[:,2],z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	gui(surf)
end

#=function plotHeatmap(x,y,z,xlab="",ylab="",name="",cmap= :inferno,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
	#surf = surface(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
	#surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	#heatmap!
	gui(surf)
end=#

function pyplotHeatmap(x,y,z,xlab="",ylab="",name="",cmap= :nipy_spectral,save=false, path="")
	fig, ax = PyPlot.subplots();
	#=
        #dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*AR
	#surf = PyPlot.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	show(size(z))
	show(size(x))
	show(size(y))
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, cmap = cmap, 
        #                      norm=matplotlib[:colors][:SymLogNorm](1e-2),
        #                     size=(width,height), margin=5mm, zscale=:log10)
        =#
        min=10^-20
        vmin = maximum([min,minimum(z)])
        surf = ax[:imshow](z.+vmin, cmap=cmap, norm=matplotlib[:colors][:LogNorm](vmin=vmin, vmax=10^-2), extent= [minimum(x), maximum(x), minimum(y), maximum(y)])
        ax[:set_xlabel](xlab)
        ax[:set_ylabel](ylab)
        PyPlot.colorbar(surf, label=name)
        if(save)
            SaveFigure(fig,path,name)
            #fig.savefig(path*name*".png")
            #close(fig)
        end
        #heatmap!
        #xlabel
	#gui(surf)
        return surf
end

function plotHeatmap(x,y,z,xlab="",ylab="",name="",cmap= :inferno,save=false, AR=0.7)
        Plots.pyplot()
	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*AR
	#surf = PyPlot.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	show(size(z))
	show(size(x))
	show(size(y))
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, cmap = cmap, 
        #                      norm=matplotlib[:colors][:SymLogNorm](1e-2),
        #                     size=(width,height), margin=5mm, zscale=:log10)
        surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height), margin=5mm, zscale=:log10)

	#surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	#heatmap!
        #xlabel
	gui(surf)
        return surf
end

function plotSurf(x,y,z,xlab="",ylab="",name="",cmap= :inferno,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	#surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
	surf = surface(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
	#surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	#heatmap!
	gui(surf)
end

#function nameGen

function Sweep2DSurf(f::Function, gridToArgs::Function, xvals::Vector, yvals::Vector, xlab="",ylab="",zlab="")
    nx = size(xvals)[1]; ny = size(yvals)[1]
    zgrid = zeros(Float64,ny,nx)
    #BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
    iter = ProgressBar(1:nx)
    for ix in iter
    #Threads.@threads for ix = 1:nx
        for iy = 1:ny
            x = xvals[ix]; y = yvals[iy];
            zgrid[iy,ix] = f(gridToArgs(x,y))
        end
    end
    plotHeatmap(xvals,yvals,zgrid,xlab,ylab,zlab,:rainbow)
end


function Sweep1DSurf(f::Function, gridToArgs::Function, xvals::Vector, yvals::Vector, xlab="",ylab="",zlab="",xscale::Float64=1,surf::Bool=true)
    zmat = zeros(size(yvals)[1],size(xvals)[1])
    stdout_global = stdout; stderr_global=stderr;
    Threads.@threads for i in eachindex(xvals)
        x = xvals[i]
        redirect_stdio(stdout=stdout_global, stderr=stderr_global)
        println("Sweeping y = $yvals at x = $x...")
        args = gridToArgs(x)
        redirect_stdio(stdout=p.path*"output.log", stderr="errors.log")
	out = f(args)
        zmat[:,i] = deepcopy(out)
        GC.gc()
        #zmat[:,i] = deepcopy(√(x)*yvals.^2)
    end
    if(surf)
        f = plotHeatmap(xvals,yvals,zmat,xlab,ylab)
    else
        f = plot1D(xscale*xvals,[zmat[i,:] for i in eachindex(yvals)],xlab, zlab,["$ylab = $y" for y in yvals],:log10)
    end
    return (xscale*xvals,reduce(hcat,[zmat[i,:] for i in eachindex(yvals)]), f)
end


function SavePlots(fig,path,name="",type=".png")
        Plots.savefig(fig, path*"/"*name*type)
        Plots.close(fig)
end


function SaveFigure(fig,path,name="",type=".png")
	fig.savefig(path*"/"*name*type)
	PyPlot.close(fig)
        GC.gc()
end

function plot3Dvectors(Rvals,fvals,cvals,xlab="x position (nm)",ylab="y position (nm)", zlab="z position (nm)", name="",cmap="bwr",scale=(1/nm),zscale=1,save=false)
	#dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	#C = 500
	#width = C
	#height = C*dy/dx
	fig = figure(figsize=(7,5))
	ax = fig.add_subplot(projection="3d");
	##fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	#zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=z);
        #maxnorm = maximum(norm.(fvals))
        #fvals = [val/maxnorm for val in fvals]
        #fvals = [val/maxnorm for val in fvals]
        R = [Tuple(scale.*reduce(vcat,Rvals')[:,i]) for i = 1:3]
        Δi = scale*[maximum(R[i]) - minimum(R[i]) for i=1:3]
        s = 15*[1,1,1]
        f = [s[i]*reduce(vcat,fvals')[:,i] for i = 1:3]
        cvals = f[2]

        #minc = minimum(cvals); maxc = maximum(cvals); cvals = (maxc-minc)^-1*(cvals.-minc)
        #cnorm = PyPlot.cm.Normalize(vmin=min.(cvals), vmax=max.(cvals))
        #fcmap = PyPlot.cm.get_cmap(cmap)
        #fcmap = PyPlot.cm.ScalarMappable(norm=cnorm, cmap=cmap)
        #vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],pivot="middle",normalize=true)
        vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],pivot="middle",normalize=true)
        #vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],pivot="middle",normalize=true,colors=fcmap(cvals))
        #vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],pivot="middle",normalize=true,colors=fcmap(cvals))
        #vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],pivot="middle",normalize=true,colors=fcmap(cvals))
        #vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3],colors=fcmap(cvals),normalize=true)
	#vectorfield = ax.quiver(R[1],R[2],R[3],f[1],f[2],f[3], cmap=cmap, c=cvals)
	#zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=Tuple([zscale*zi for zi in z]), cmap=cmap);
        xmin = minimum(R[1]); xmax = maximum(R[1])
        ymin = minimum(R[2]); ymax = maximum(R[2])
        max = maximum([ymax,xmax]); min = maximum([ymin,xmin])
        PyPlot.xlim(min,max)
        PyPlot.ylim(-max,max)
        PyPlot.xlabel(xlab);
        PyPlot.ylabel(ylab);
        PyPlot.zlabel(zlab);
	#fig.colorbar(zplot, ax=ax);
	#PyPlot.colorbar(vectorfield, label=name);
        pygui()
		#pygui(vectorfield)
	#gcf()
	#show(surf);
end

function plotScatter(R,z,xlab="",ylab="",name="",cmap="inferno",xyscale=(1/nm),zscale=1,save=false)
	dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	C = 500
	width = C
	height = C*dy/dx
	fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	#zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=z);
        max = maximum(abs.([maximum(z),minimum(z)]))
	zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=Tuple([zscale*zi for zi in z]), cmap=cmap, vmin = -max, vmax = max);
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


function plotBands(klist, nk, E, max::String="auto")
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
        if(max=="auto")
            maxE = maximum(E)
            minE = minimum(E)
            ylim(minE,maxE)
        else
            maxE = Float64(max)
            minE = -2.0
            ylim(minE,maxE)
        end
        #maxE = Emax
	#minE = -Emax
	ylabel("E - Ef (eV)")
	for iE = 1:nE
		Evals = collect(E[:,iE])
		#plot(indices,Evals)
		PyPlot.scatter(indices,Evals,s=1)
		#display(plot!(indices,Evals))
	end
	gcf()
end


for n in names(@__MODULE__; all=true)
   if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
       @eval export $n
   end
end


end
