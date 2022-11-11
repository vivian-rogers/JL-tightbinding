

module DOS
using UsefulFunctions
#using Arpack
using LinearAlgebra
using Arpack
using GMaterialParameters
using Logging
using PlotlyJS
using Distributed
using ProgressBars
using StatsBase
using Operators
#using GLMakie

export energySurface, kslice, getDOS, getLDOS, eigSurface


function genBZmesh(B::Matrix,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    X1 = B*[1/2;0;0];
    X2 = B*[0;1/2;0];
    X3 = B*[0;0;1/2];
    kx =  LinRange(-X1[1],X1[1],2*nx+1)
    ky =  LinRange(-X2[2],X2[2],2*ny+1)
    kz =  LinRange(-X3[3],X3[3],2*nz+1)
    return kx, ky, kz
end


function kDOSgen(p::NamedTuple,H::Function,E::Float64)
    Γ = [1;1;1]*0.0
    testH = H(Γ)
    n = size(testH)[1]
    function kDOS(k::Vector{Float64})
        # approximate the GF
        #=
        λs, eigstates = eigs(H(k).-(E+im*10^-4)*I(n),nev=1,which=:SM)
        display(λ)
        Gᴿ = sum(map(λ-> 1/(im*p.η - λ), λs))
        return real(tr(Gʳ))
        =#
        # estimate a bound for the condition number
        
        #diagonal = diag(H(k)-(E)*I(n))
        #return maximum(abs.(diagonal))/minimum(abs.(diagonal))
        
        # Plot the real Gʳ for the 0-crossing
        #return real(tr(grInv((E+p.η)*I(n) .- H(k))))
        
        # Plot the DOS
        return imag(tr(grInv((E+p.η)*I(n) .- H(k))))
    end
    return kDOS
end

e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1]


function kslice(p::NamedTuple,H::Function, E::Float64, kfixed::String="x", nk1::Int=200, nk2::Int=200, k3::Float64=0)
        DOS = zeros(2*nk1+1,2*nk2+1)
        nx = ny = nz = 1
        ik = 0
        if(kfixed=="x")
            ik = 1;
            ny = nk1; nz = nk2
        elseif(kfixed=="y")
            ik = 2;
            nx = nk1; nz = nk2
        else
            ik = 3;
            nx = nk1; ny = nk2
        end
        X₃ = zeros(3); X₃[ik] = 1/2; X₃ = p.B*X₃; 
        ks = genBZmesh(p.B,nx,ny,nz)
        iterdims = setdiff([1,2,3],[ik])
        k1 = collect(ks[minimum(iterdims)]); k2 = collect(ks[maximum(iterdims)]);
        kDOS = kDOSgen(p,H,E)
        println("Imaging $kfixed slice DOS($(E) eV) for $((nk1*2+1)*(nk2*2+1)) grid")
        iter = ProgressBar(1:(2*nk1+1))
        #for ik1 in iter
        #for ik1 = 1:(2*nk1+1)
        BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
        Threads.@threads for ik1 in iter
            for ik2 = 1:(2*nk2+1)
                k = zeros(3); k[ik] = X₃[ik]; 
                k[minimum(iterdims)] = k1[ik1]; k[maximum(iterdims)] = k2[ik2]; 
                DOS[ik1,ik2] = kDOS(k)
            end
        end
        f = plot(heatmap(
             x = k1,
             y = k2,
             z = log10.(DOS)
        ))
         display(f)
         #plotHeatmap(k1,k2,DOS',"k₂","k₃","DOS(ky,kz)",:)
       return f 
end


function eigSurface(p::NamedTuple,H::Function, Q::Matrix, neigs::Int=4, kfixed::String="x", nk1::Int=200, nk2::Int=200, k3::Float64=0, cones::Bool=true)
        eigsheets = zeros(neigs,2*nk1+1,2*nk2+1)
        projectionsheets = zeros(neigs,2*nk1+1,2*nk2+1)
        conesheets = zeros(neigs,2*nk1+1,2*nk2+1,6)
        nx = ny = nz = 1
        ik = ik₃ = ik₂ = 0
        if(kfixed=="x")
            ik = 1;
            ik₂ = 2; ik₃ = 3
            ny = nk1; nz = nk2
        elseif(kfixed=="y")
            ik = 2;
            ik₂ = 1; ik₃ = 3;
            nx = nk1; nz = nk2
        else
            ik = 3;
            ik₂ = 1; ik₃ = 2;
            nx = nk1; ny = nk2
        end
        X₃ = zeros(3); X₃[ik] = 1/2; X₃ = p.B*X₃; 
        ks = genBZmesh(p.B,nx,ny,nz)
        iterdims = setdiff([1,2,3],[ik])
        k1 = collect(ks[minimum(iterdims)]); k2 = collect(ks[maximum(iterdims)]);
        println("Imaging $kfixed slice spectrum for $neigs eigenvalues, $((nk1*2+1)*(nk2*2+1)) grid")
        iter = ProgressBar(1:(2*nk1+1))
        #for ik1 in iter
        #for ik1 = 1:(2*nk1+1)
        #BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
        n = size(H([0.0;0.0;0.0]))[1]
        nospinN = Int(n/2)
        #Sx = I(nospinN)⊗σ₁; Sy = I(nospinN)⊗σ₂; Sz = I(nospinN)⊗σ₃
        for ik1 in iter
        #Threads.@threads for ik1 in iter
            for ik2 = 1:(2*nk2+1)
                k = zeros(3); k[ik] = k3*X₃[ik]; 
                k[minimum(iterdims)] = k1[ik1]; k[maximum(iterdims)] = k2[ik2]; 
                Energies, Eigenstates = eigs(H(k), nev = neigs, which=:SM)
                sortE = sortperm(real.(Energies))
                Energies = Energies[sortE]
                Eigenstates = Eigenstates[:,sortE] 
                #Energies = real.(Energies(sortE))
                eigsheets[:,ik1,ik2] = real.(Energies)
                for iE = 1:neigs
                    eigstate = Eigenstates[:,iE]
                    projectionsheets[iE,ik1,ik2] = real(eigstate'*Q*eigstate)
                    if(cones)
                        spinperm = [ik₂;ik₃; ik]
                        for ax = 1:3
                            Sᵢ = I(nospinN)⊗σ[spinperm[ax]] 
                            conesheets[iE,ik1,ik2,ax] = real(eigstate'*Sᵢ*eigstate)
                        end
                        conesheets[iE,ik1,ik2,4] = k1[ik1]
                        conesheets[iE,ik1,ik2,5] = k2[ik2]
                        conesheets[iE,ik1,ik2,6] = real(Energies[iE])
                    end
                end
            end
        end
        # to get the effective edge of conduction, valence bands
        minmax = minimum(map(maximum∘(f(i1,i2) = eigsheets[:,i1,i2]),collect(1:(2*nk1+1)),collect(1:(2*nk2+1))))
        maxmin = maximum(map(minimum∘(f(i1,i2) = eigsheets[:,i1,i2]),collect(1:(2*nk1+1)),collect(1:(2*nk2+1))))
        maxE = minmax*1.2+1.1; minE = maxmin*1.2-1.1;
        layout = Layout(scene = attr(
                                     xaxis_title = "k₁ (1/m)",
                                     yaxis_title = "k₂ (1/m)",
                                     zaxis_title = "Energy (eV)",
                                     zaxis=attr(range=[minE,maxE])
                                     ),
                        title="σₓ-projected band structure for fixed k$kfixed = $(round(k3,sigdigits=2))×X ",)
        if(cones)
            surfaces = []
            surfaces = [surface(x=k1, y=k2, z=eigsheets[iE,:,:], surfacecolor=projectionsheets[iE,:,:], 
                            colorscale=colors.RdBu_3, opacity=0.5, showscale = (iE==1), cmin = minimum(projectionsheets), cmax = maximum(projectionsheets))
                            for iE = 1:neigs] 
            x = vec(conesheets[1,:,:,4]); y = vec(conesheets[1,:,:,5]); z = vec(conesheets[1,:,:,6]);
            u = vec(conesheets[1,:,:,1]); v = vec(conesheets[1,:,:,2]); w = vec(conesheets[1,:,:,3]);
            dE = maxE - minE; dk1 = maximum(k1) - minimum(k1); dk2 = maximum(k2) - minimum(k2)
            S = sx = sy = sz = 1
            
            #S = 10^-3*dk1*dk2*0.2;
            #sx = S/(dE*dk2); sy = S/(dE*dk1); sz = S/(dk1*dk2); 
            coneplots = [cone(x = vec(conesheets[iE,:,:,4]), y = vec(conesheets[iE,:,:,5]), z = vec(conesheets[iE,:,:,6]),
                                u = sx*vec(conesheets[iE,:,:,1]), v = sy*vec(conesheets[iE,:,:,2]), w = sz*vec(conesheets[iE,:,:,3]),
                                colorscale = [[0, "rgb(120,120,120)"],[1,"rgb(120,120,120)"]], sizemode="scaled", 
                                showscale=false, hoverinfo="u+v+w+name", sizeref=0.1) for iE = 1:neigs]
            #iE = 1
            #S = 10^6
            #xt = k1; yt = k2; zt = vec(eigsheets[iE,:,:]); 
            #f2 = plot(cone(x=xt,y=yt,z=zt, u=ut, 
                          #v=vt, w=wt, sizemode="absolute"), layout)
            #v=vt, w=wt, sizemode="scaled", sizeref=10^6))
            #display(f2)
            #f = plot(cone(x=xt,y=yt,z=zt, u=ut, 
                          #v=vt, w=wt, sizemode="absolute"), layout)
            #v=vt, w=wt, sizemode="scaled", sizeref=0.25, colorscale=[[0, "rgb(120,120,120)"], [1, "rgb(120,120,120)"]]))
            #f = plot(coneplots)
            f = plot(vcat(surfaces,coneplots),layout)
            display(f)
           return f 
       end
        surfaces = [surface(x=k1, y=k2, z=eigsheets[iE,:,:], surfacecolor=projectionsheets[iE,:,:], 
                            colorscale=colors.RdBu_3, opacity=0.5, showscale = (iE==1), cmin = minimum(projectionsheets), cmax = maximum(projectionsheets))
                            for iE = 1:neigs] 
        f = plot(surfaces,layout)
        display(f)
       return f 
end



function energySurface(p::NamedTuple,H::Function, E::Float64, nx::Int=10, ny::Int=10, nz::Int=10)
        kx, ky, kz = genBZmesh(p.B,nx,ny,nz)
        kDOS = kDOSgen(p,H,E)
        println("Evaluating DOS over brillouin zone, may take a while")
        # for plotlyjs
        S = 1
        X1 = S*p.B*[1/2;0;0];
        X2 = S*p.B*[0;1/2;0];
        X3 = S*p.B*[0;0;1/2];
        k1 =  LinRange(-X1[1],X1[1],2*nx+1)
        k2 =  LinRange(-X2[2],X2[2],2*ny+1)
        k3 =  LinRange(-X3[3],X3[3],2*nz+1)
        kx, ky, kz = mgrid(k1,k2,k3)
        # for makie
        #kDOSmesh = [kDOS(kx*e1+ky*e2+kz*e3) for kx = kx, ky = ky, kz = kz]
        BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
        #f(kx::Float64,ky::Float64,kz::Float64) = tanh(kDOS(e1*kx + e2*ky + e3*kz))
        kvals = CartProd([collect(k3),collect(k2),collect(k1)])
        #display(kvals)
        vals = zeros(size(kvals)[1])
        iter = ProgressBar(eachindex(kvals))
        Threads.@threads for ik in iter
            k = kvals[ik][:]
            #display(k)
            vals[ik] = kDOS(k)
        end
        vals = log.(vals)
        vals = (1/maximum(vals))*(vals .- minimum(vals)) #normalize to maximum of 3 for tanh
        vals = tanh.(vals)
        f1 = plot(histogram(x=vals))
        display(f1)
        kDOSmesh = reshape(vals,(nx*2+1,ny*2+1,nz*2+1))
        # may not be optimal but it works
        #kDOSmesh = map(f, kx,ky,kz)
        #kDOSmesh = pmap(k1,k2,k3 -> kDOS(k1*e1+k2*e2+k3*e3), kx
        println("Done! Plotting energy surface")
        ϵ = 10^-7; 
        #max = ϵ; min = -ϵ;
        max = percentile(vals,99.8); min = percentile(vals,98.5)
        f = plot(isosurface(
            x=kx[:],
            y=ky[:],
            z=kz[:],
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
        ))
        #f = volume(kDOSmesh,algorithm = :iso, isorange = 0.02, isovalue=0.5)
        #show(f)
        #display(kx)
        display(f)
        #gui()
        return f
end

#=function getLDOS(kgrid, weights, Hofk, λ, Emax=4, Emin=-4, Ebins=100) #takes in array of kpt strings, number of interpolation pts, H(k) function
	testH = Hofk([0;0])
	nk = size(kgrid)[1]
	#initialize the band array
	#λ_test, evecs_test = eig(testH, 1E-12)
	N = 2*λ^2
	maxiter = 400
	nE = 12*Int(floor(log2(size(testH)[1])))
	nEig = size(testH)[1]
	if(nE < size(testH)[1])
		println("Heads up! $nE / $nEig eigvls are being calculated")
	end
	λ_test, evecs_test = eigs(testH,nev=nE, maxiter=maxiter)
	ndim = size(evecs_test)[2]
	LDOS = zeros(Ebins,2*λ+1)
	DOS = zeros(Ebins)
	Evals = LinRange(Emin,Emax,Ebins)
	function idos(E)
		if(E < Emax && E>Emin)
			return Int(round((E-Emin)*Ebins/(Emax-Emin)))
		else
			return false
		end
	end
	Logging.disable_logging(Logging.Warn)
	for ik in 1:nk
		k = kgrid[ik]
		w = weights[ik]
		H = Hofk(k)
		Eofk, Estatek = eigs(H,nev=nE, which=:SM, sigma=10^-5, maxiter=maxiter)
		for iE in 1:nE
			Eval = real(Eofk[iE]) # particular energy value
			Estate = Estatek[:,iE]
			index = idos(Eval)
			if(index != false)
				# add to real DOS
				DOS[index] += w
				#take a diagonal path!
				for i = 1:λ
					i1 = g2i(i,i,λ)
					overlap = abs(Estate[i1])^2
					LDOS[index,2*i-1] += w*overlap
				end
				for i = 1:λ
					i2 = g2i(i+1,i,λ)
					overlap = abs(Estate[i2])^2
					LDOS[index,2*i] += w*overlap
				end
			end
		end
	end
	return LDOS, DOS, Evals
end
function getDOS(kgrid, weights, Hofk, Emax=1, Emin=-1, Ebins=100) #takes in array of kpt strings, number of interpolation pts, H(k) function
	testH = Hofk([0;0])
	nk = size(kgrid)[1]
	#initialize the band array
	#λ_test, evecs_test = eig(testH, 1E-12)
	maxiter = 400
	nE = 4*Int(floor(log2(size(testH)[1])))
	nEig = size(testH)[1]
	if(nE < size(testH)[1])
		println("Heads up! $nE / $nEig eigvls are being calculated")
	end
	λ_test, evecs_test = eigs(testH,nev=nE, maxiter=maxiter)
	ndim = size(evecs_test)[2]
	DOS = zeros(Ebins)
	Evals = LinRange(Emin,Emax,Ebins)
	function idos(E)
		if(E <= Emax && E>=Emin)
			return Int(round((E-Emin)*Ebins/(Emax-Emin)))
		else
			return false
		end
	end
	
	# I do this so that I am not yelled at by the eigensolver for shifting the spectrum
	Logging.disable_logging(Logging.Warn)
	for ik in 1:nk
		k = kgrid[ik]
		w = weights[ik]
		H = Hofk(k)
		Eofk, Estatek = eigs(H,nev=nE, which=:SM, sigma=10^-5, maxiter=maxiter)
		for iE in 1:nE
			Eval = real(Eofk[iE]) # particular energy value
			index = idos(Eval)
			if(index != false)
				DOS[index] += w
			end
		end
	end
	return DOS, Evals
end=#

end
