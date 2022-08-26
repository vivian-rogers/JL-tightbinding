module MakiePlots

using Makie
using GLMakie
#using GeometryTypes
using LinearAlgebra
using Constants
using ColorSchemes
using UsefulFunctions

export render2TDevice

function slab(xlims::Vector{Float64}, ylims::Vector{Float64}, zlims::Vector{Float64})
    dx = xlims[2] - xlims[1]; dy = ylims[2] - ylims[1]; dz = zlims[2] - zlims[1];
    return FRect3D( Point3f0(avg(xlims),avg(ylims),avg(zlims)), Point3f0(dx,dy,dz))
end


# make plotting easier on the Spirit
nm = 1
convnm = 10^9

function plotslab!(slab, color::Symbol, alpha::Float64=0.5)
    mesh!(slab, color=(color,alpha), transparency=true)
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

function render2TDevice(p::NamedTuple,Rvals::Vector{Vector{Float64}}, fvals::Vector{Vector{Float64}}, cmap::Function=norm(), margin::Float64=2*nm, arrowScaling::Float64=1.0)
        # make the arrow plot
        println("Rendering 2 Terminal device with magnetization field")
        Rvals = Rvals*convnm
        points = Point3f0.(Rvals)
        arrowvec = 0.1*Vec3f0.(fvals)
        colors = cmap.(fvals)
        display(points)
        display(arrowvec)
        fig = arrows(points, arrowvec, fxaa=true, color=colors, linewidth=0.1, arrowsize=Vec3f0(0.3,0.3,0.4), axis=(type=Axis3,), align = :center)
        # start setting up the rendering for the boxes
        dx = p.a*convnm
        xmin = -dx/2; xmax = convnm*p.SLa₁[1]-dx/2
        ymin = minimum([Rval[2] for Rval in Rvals])-margin
        ymax = maximum([Rval[2] for Rval in Rvals])+margin
        zmin = minimum([Rval[3] for Rval in Rvals])-margin
        zmax = maximum([Rval[3] for Rval in Rvals])+margin
        totlim = 10*nm
        Lcontact = slab([-totlim+xmin,xmin],[ymin,ymax],[zmin,zmax])
        device = slab([xmin,xmax],[ymin,ymax],[zmin,zmax])
        Rcontact = slab([xmax,totlim+xmax],[ymin,ymax],[zmin,zmax])
        #wireframe!([Lcontact,device,Rcontact])
        alpha = 0.5
        plotslab!(Lcontact,:grey,alpha)
        plotslab!(device,:green,alpha)
        plotslab!(Rcontact,:grey,alpha)
        display(fig)
end

end
