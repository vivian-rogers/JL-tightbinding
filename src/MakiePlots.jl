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
using ColorSchemes
using UsefulFunctions

export render2TDevice

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
        fig = Figure(figure_padding=50,resolution=(1400,800))
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
