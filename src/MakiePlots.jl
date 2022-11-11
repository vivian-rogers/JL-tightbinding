module MakiePlots
using Makie
#using RPRMakie
using GLMakie
using GLMakie.GLFW
using GLMakie: to_native
#using GeometryTypes
using LinearAlgebra
using Constants
#using Colors
using Distributed
using StatsBase
using ColorSchemes
using UsefulFunctions
using ProgressBars
using PlotStuff

export render2TDevice, mixedDOS

#=RPRMakie.activate!(;
    iterations=200,
    resource=RPR.RPR_CREATION_FLAGS_ENABLE_GPU0,
    plugin=RPR.Tahoe)
=#

function slab(xlims::Vector{Float64}, ylims::Vector{Float64}, zlims::Vector{Float64})
    dx = xlims[2] - xlims[1]; dy = ylims[2] - ylims[1]; dz = zlims[2] - zlims[1];
    #return FRect3D( Point3f0(xlims[1],ylims[1],zlims[1]), Point3f0(xlims[2],ylims[2],zlims[2]))
    return Makie.FRect3D( Point3f(xlims[1],ylims[1],zlims[1]), Point3f(dx,dy,dz))
    #return FRect3D( Point3f0(avg(xlims),avg(ylims),avg(zlims)), Point3f0(dx,dy,dz))
end


# make plotting easier on the Spirit
nm = 1

function plotslab!(scene, slab, color::Symbol, alpha::Float64=0.5)
    mesh!(scene, slab, color=(color,alpha), transparency=true)
    #=points = Makie.decompose(Makie.Point3f0, slab)
    faces = Makie.decompose(GLTriangle, slab)
    mesh!(
        points, faces,
        color = LinRange(1.0,1.0, 8),
        colormap = [
            (color, alpha),
        ],
        transparency = true
    )=#
end

function render2TDevice(p::NamedTuple,Rvals::Vector{Vector{Float64}}, fvals::Vector{Vector{Float64}}, cmap::Function=norm(), margin::Float64=5.0, arrowScaling::Float64=1.0)
        # make the arrow plot
        convnm = 10^9
        println("Rendering 2 Terminal device with magnetization field")
        dx = p.a*convnm
        xmin = -dx/2; xmax = convnm*p.SLa₁[1]-dx/2
        ymin = minimum([Rval[2] for Rval in Rvals])-margin
        ymax = maximum([Rval[2] for Rval in Rvals])+margin
        zmin = minimum([Rval[3] for Rval in Rvals])-margin
        zmax = maximum([Rval[3] for Rval in Rvals])+margin
        totlim = 5*nm+margin
        ϵ = 0.0000001*nm
        
        Rvals = Rvals*convnm
        colors = cmap.(fvals)
        cmax = maximum(abs.(colors));
        Rvals = push!(Rvals,[-100.0;0.0;0.0]); 
        Rvals = push!(Rvals,[-100.0;0.0;0.0]); 
        fvals = push!(fvals,[0;cmax;0]); 
        fvals = push!(fvals,[0;-cmax;0]); 
        colors = push!(colors,cmax)
        colors = push!(colors,-cmax)
        points = Point3f.(Rvals)
        lengths=norm.(fvals)
        maxlength=maximum(lengths)
        arrowvec = 0.5*(margin/maxlength)*Vec3f.(fvals)
        #display(points)
        #display(arrowvec)
        # Testing with arrow and a box
        #fig = arrows([Point3f0(0)],[Vec3f0(1.5)],linewidth=0.1, axis=(type=Axis3,))
        #shape = FRect((0, 0, 0), (2, 2, 2))
        #color = RGBAf0(.2, .8, .2, .4)
        #mesh!(shape; color, transparency=true)
        
        
        #fig = arrows(points, 10^-12*arrowvec, fxaa=true, color=colors, linewidth=0.1, arrowsize=Vec3f0(0.3,0.3,0.4), axis=(type=Axis3,), align = :center)
        # start setting up the rendering for the boxes
        Lcontact = slab([-totlim+xmin,xmin],[ymin,ymax],[zmin,zmax])
        device = slab([xmin+ϵ,xmax-ϵ],[ymin,ymax],[zmin,zmax])
        Rcontact = slab([xmax,totlim+xmax],[ymin,ymax],[zmin,zmax])
        #wireframe!([Lcontact,device,Rcontact])
        lim = slab([-totlim+xmin,totlim+xmax],[ymin,ymax],[zmin,zmax])
        fig = Figure(figure_padding=50,resolution=(1000,500))
        fontsize_theme = Theme(fontsize = 25)
        set_theme!(fontsize_theme)
        aspect = ((xmax+totlim - xmin-totlim)/(ymax-ymin),1,1)
        #smalllabel=20;
        scene = Axis3(fig[1,1],aspect = aspect, xlabel = "X Position (nm)", ylabel = "Periodic in Y",
                      zlabel = "Periodic in Z", title = "DW start position = $(round(p.startDWs*convnm,sigdigits=3)) nm")
        
        alpha = 0.3
        scene.xticks = 0.0:20*dx:(xmax+dx)
        scene.yticklabelsvisible = false
        scene.yticksvisible = false
        scene.zticklabelsvisible = false
        scene.zticksvisible = false
        plotslab!(scene, Lcontact,:grey,alpha)
        plotslab!(scene, device,:green,alpha)
        plotslab!(scene, Rcontact,:grey,alpha)
        cmap = :coolwarm; 
        #cmap = colormap("RdBu", 100; mid=0.0)
        arrows!(scene, points, arrowvec, fxaa=true, color=colors, linewidth=0.3, arrowsize=Vec3f(0.6,0.6,0.8),
                colormap=cmap,  align = :origin)
        Colorbar(fig[1, 2], limits = (-maximum(abs.(colors)),maximum(abs.(colors))), colormap = cmap, height = Relative(0.5), label = "β⋅ŷ (eV)") 
        limits!(scene,lim)
        #scene.xlabel("X Position (nm)") 
        #scene.ylabel("Periodic in Y")
        #scene.zlabel("Periodic in Z") 
        #scene.title("DW start position = $(round(p.startDWs*convnm,sigdigits=3)) nm")
        #arrows!(points, arrowvec, fxaa=true, color=colors, linewidth=0.05, arrowsize=Vec3f0(0.1,0.1,0.15), axis=(type=Axis3,), limits=lim, align = :center)
        #limits!(fig,lim)
        if(p.savedata)
            save(p.path * "2TDeviceRender.png",fig)
            #glfw_window = to_native(display(scene))
            #GLFW.SetWindowShouldClose(glfw_window,true)
        else
            display(fig)
        end
end

function mixedDOS(p::NamedTuple,DOS::Function, ny::Int=10, nz::Int=10)
        #kx, ky, kz = genBZmesh(p.B,nx,ny,nz)
        #kDOS = kDOSgen(p,H,E)
        #DOS_Γ = DOS([0;0;0])
        
        println("Evaluating DOS over brillouin zone, may take a while")
        # for plotlyjs
        S = 0.15; 
        X2 = S*p.B*[0;1/2;0];
        X3 = S*p.B*[0;0;1/2];
        k2 =  LinRange(-X2[2],X2[2],2*ny+1)
        k3 =  LinRange(-X3[3],X3[3],2*nz+1)
        ixs =  LinRange(1,p.nx,p.nx)
        #ky, kz, xs = mgrid(k2,k3,ixs)
        # for makie
        #kDOSmesh = [kDOS(kx*e1+ky*e2+kz*e3) for kx = kx, ky = ky, kz = kz]
        BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
        #f(kx::Float64,ky::Float64,kz::Float64) = tanh(kDOS(e1*kx + e2*ky + e3*kz))
        #kvals = CartProd([collect(k2),collect(k3)])
        #display(kvals)
        #vals = zeros(p.nx*ny*nz)
        #vals = zeros(size(kvals)[1])
        kDOSmeshL = zeros(2*ny+1,2*nz+1,p.nx)
        kDOSmeshR = zeros(2*ny+1,2*nz+1,p.nx)
        #=iter = ProgressBar(eachindex(k2))
        Threads.@threads for ik2 in iter=#
        println("Running fermi surface over $(2*ny+1) x $(2*nz+1) k-grid")
        for ik2 in eachindex(k2)
            ky = k2[ik2]
            for ik3 in eachindex(k3)
                kz = k3[ik3]
                k = [0,ky,kz]
                D = DOS(k)
                kDOSmeshL[ik2,ik3,:] = D[1]
                kDOSmeshR[ik2,ik3,:] = D[2]
            end
        end
        #DOS = zeros(p.nx*(2*ny+1)*(2*nz+1))
        #for iDOS = eachindex(kDOSmesh)
        #    DOS[iDOS] = kDOSmesh[iDOS]
        #end
        #kDOSmesh = log10.(kDOSmesh)
        fig = Figure(figure_padding=40,resolution=(1000,500))
        cmap = :inferno
        #n = 101
        #g(x) = x^2
        #alphas = [g(x) for x in range(0, 1, length = n)]
        #cmap_alpha = resample_cmap(cmap, n; alpha = alphas)
        
        xs = LinRange(-S,S,2*ny+1); ys = LinRange(-S,S,2*nz+1); zs = collect(1:p.nx)
        points3d = [Point3f(ix, iy, iz) for ix in xs, iy in ys, iz in zs];
        #points = []; vals = []
        #cutoff = 3; points = vec(points3d); DOS = vec(kDOSmesh);
        #for DOS in vec(kDOSmesh)
        #    if(DOS > cutoff)
        #        points
        #scatter!(scene, points[DOS .> cutoff]; color = DOS[DOS .> cutoff], colormap=cmap)
        pc(percent) = percentile(vec(kDOSmeshR),percent)
        println("99.9% = $(pc(99.9)), 99.5% = $(pc(99.5)), 99.0% = $(pc(99.0)), 97% = $(pc(97.0))")
        #contour!(scene, LinRange(-S,S,2*ny+1),LinRange(-S,S,2*nz+1),collect(1:p.nx), 
        #         kDOSmesh, alpha=0.2, levels = LinRange(pc(99.0),pc(99.9),5), colormap=:viridis)
        #fig2 = Figure(); ax2 = Axis3(fig[1,1]); 
        #hist(fig[1,2], vec(kDOSmesh), bins = 100);
        #fh = hist(vec(kDOSmesh), bins = 100)
        #display(fh)
        scene = Axis3(fig[1,1], aspect = (1,1,1), title="PDOS(H₀+Σₗ+Σᵣ)",
                      xlabel = "ky (2π/a)", ylabel = "kz (2π/a)", zlabel = "X Position (nm)")
        #scene2 = Axis3(fig[1,2], aspect = (1,1,1), title="¼⟨1-γ⁵⟩⟨1+γ⁵⟩",
        scene2 = Axis3(fig[1,2], aspect = (1,1,1), title="⟨L⟩+⟨R⟩",
                      xlabel = "ky (2π/a)", ylabel = "kz (2π/a)", zlabel = "X Position (nm)")
        volume!(scene, LinRange(-S,S,2*ny+1),LinRange(-S,S,2*nz+1),collect(1:p.nx), 
                kDOSmeshL, algorithm=:iso, alpha=0.7, isorange = 0.1, isovalue = 0.3, colormap=:redsblues)
        volume!(scene, LinRange(-S,S,2*ny+1),LinRange(-S,S,2*nz+1),collect(1:p.nx), 
                kDOSmeshR, algorithm=:iso, alpha=0.7, isorange = 0.1, isovalue = 0.3, colormap=:bluesreds)
        both = kDOSmeshR.+kDOSmeshL
        volume!(scene2, LinRange(-S,S,2*ny+1),LinRange(-S,S,2*nz+1),collect(1:p.nx), 
                both, algorithm=:iso, alpha=0.7, isorange = percentile(vec(both),98.0)/2, isovalue = percentile(vec(both),98.0), colormap=:spring)
        #scale!(fig,1.0,1.0,1.0)
        if(p.savedata)
            save(p.path * "MixedDOS.png",fig)
            #glfw_window = to_native(display(scene))
            #GLFW.SetWindowShouldClose(glfw_window,true)
        else
            display(fig);
        end
        #vals = vec(kDOSmesh) #add back in Log here
        #vals = (1/maximum(vals))*(vals .- minimum(vals)) #normalize to maximum of 3 for tanh
        #vals = tanh.(vals)
        #f1 = plot(histogram(x=vals))
        #display(f1)
        #kDOSmesh = reshape(vals,(p.nx,ny*2+1,nz*2+1))
        # may not be optimal but it works
        #kDOSmesh = map(f, kx,ky,kz)
        #kDOSmesh = pmap(k1,k2,k3 -> kDOS(k1*e1+k2*e2+k3*e3), kx
        #println("Done! Plotting energy surface")
        #ϵ = 10^-7; 
        #max = ϵ; min = -ϵ;
        #max = percentile(vals,99.8); min = percentile(vals,98.5)
        #=f = plot(isosurface(
            x=collect(1:p.nx),
            y=k2[:],
            z=k3[:],
            value=kDOSmesh[:],
            isomin=min,
            isomax=max,
            #=colorscale=[
                    [0, "rgb(240, 237, 74)"],
                    [1.0, "rgb(240, 237, 74)"]
            ],=#
            opacity=0.3, # needs to be small to see through all surfaces
            #opacityscale="max",
            surface_count=3, # needs to be a large number for good volume rendering
        ))=#
        #f = volume(kDOSmesh,algorithm = :iso, isorange = 0.02, isovalue=0.5)
        #show(f)
        #display(kx)
        #display(f)
        #gui()
        return fig
end
#=
function render2TDevice(p::NamedTuple,Rvals::Vector{Vector{Float64}}, fvals::Vector{Vector{Float64}}, cmap::Function=norm(), margin::Float64=5, arrowScaling::Float64=1.0)
        # make the arrow plot
        println("Rendering 2 Terminal device with magnetization field")
        Rvals = Rvals*convnm
        points = Point3f.(Rvals)
        lengths=norm.(fvals)
        maxlength=maximum(lengths)
        arrowvec = 0.7*(margin/maxlength)*Vec3f.(fvals)
        colors = cmap.(fvals)
        #display(points)
        #display(arrowvec)
        # Testing with arrow and a box
        #fig = arrows([Point3f0(0)],[Vec3f0(1.5)],linewidth=0.1, axis=(type=Axis3,))
        #shape = FRect((0, 0, 0), (2, 2, 2))
        #color = RGBAf0(.2, .8, .2, .4)
        #mesh!(shape; color, transparency=true)
        
        
        #fig = arrows(points, 10^-12*arrowvec, fxaa=true, color=colors, linewidth=0.1, arrowsize=Vec3f0(0.3,0.3,0.4), axis=(type=Axis3,), align = :center)
        # start setting up the rendering for the boxes
        dx = p.a*convnm
        xmin = -dx/2; xmax = convnm*p.SLa₁[1]-dx/2
        ymin = minimum([Rval[2] for Rval in Rvals])-margin
        ymax = maximum([Rval[2] for Rval in Rvals])+margin
        zmin = minimum([Rval[3] for Rval in Rvals])-margin
        zmax = maximum([Rval[3] for Rval in Rvals])+margin
        totlim = 5*nm+margin
        ϵ = 0.0000001*nm
        Lcontact = slab([-totlim+xmin,xmin],[ymin,ymax],[zmin,zmax])
        device = slab([xmin+ϵ,xmax-ϵ],[ymin,ymax],[zmin,zmax])
        Rcontact = slab([xmax,totlim+xmax],[ymin,ymax],[zmin,zmax])
        #wireframe!([Lcontact,device,Rcontact])
        lim = slab([-totlim+xmin,totlim+xmax],[ymin,ymax],[zmin,zmax])
        fig = Figure()
        fontsize_theme = Theme(fontsize = 35)
        set_theme!(fontsize_theme)
        scene = LScene(fig[1,1], scenekw = (limits=lim, camera=cam3d!, show_axis=true, 
                xlabel = "X Position (nm)", ylabel = "Periodic in Y", zlabel = "Periodic in Z", title = "DW start position = $(round(p.startDWs*convnm,sigdigits=3)) nm"
               ))
        alpha = 0.3
        plotslab!(scene, Lcontact,:grey,alpha)
        plotslab!(scene, device,:green,alpha)
        plotslab!(scene, Rcontact,:grey,alpha)
        arrows!(scene, points, arrowvec, fxaa=true, color=colors, linewidth=0.3, arrowsize=Vec3f(0.45,0.45,0.6), 
                xlabel = "X Position (nm)", ylabel = "Periodic in Y", zlabel = "Periodic in Z", title = "DW start position = $(round(p.startDWs*convnm,sigdigits=3)) nm",
                colormap=:coolwarm,  align = :origin)
        #arrows!(points, arrowvec, fxaa=true, color=colors, linewidth=0.05, arrowsize=Vec3f0(0.1,0.1,0.15), axis=(type=Axis3,), limits=lim, align = :center)
        #limits!(fig,lim)
        display(fig)
end
=#
end
