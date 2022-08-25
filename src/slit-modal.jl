### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 3db88667-117f-4485-b477-a6c886ce8c52
begin
    import Pkg; Pkg.activate(ENV["PWD"])
	Pkg.add("Plots")
	Pkg.add("Contour")
end

# ╔═╡ a2228ad2-61da-4953-94a5-a5585b5c3a05
begin
	using LinearAlgebra
	import Plots
	import Contour
end

# ╔═╡ 8972d564-4a82-49b7-b9bc-b97812aa7f05
implus(z) = imag(z) < 0 ? -z : z

# ╔═╡ 7f4ee055-69aa-4603-a82c-d86e616721eb
md"""
I consider a rectangular mesh of t,x,z of intervals dt,dx,dz with Nt,Nx,Nz meshes

The excitation is monochromatic, the wavelength is 1.

The period is d = 1/2. There are therefore Nx = Nt/2 meshes in one period.

The slit width ws is d/6. Therefore Ns = Nx/6

Cf. Direct computation of ZEFSs Problem with boundary conditions surfaces perpendicular to the conducting surfaces.
"""


# ╔═╡ 6a0001ef-699c-4669-8c0e-7cac66f40cc9
md"""
### 7 points stencil
"""

# ╔═╡ fb6c7048-9566-4700-9c2b-9c1957225be6
function fkz_7pts(dt, dx, dz, kxN)
    phi = (cos(dt*2pi)-1)*(dz/dt)^2 - (cos(kxN*dx*2pi)-1)*(dz/dx)^2
    implus(acos(1.0 + 0.0im + phi)/(dz*2pi))
end

# ╔═╡ 9cdbc786-3f96-4dce-bc9d-d3374ec7ab91
fkz_7pts(delta, kxN) = fkz_7pts(delta, delta, delta, kxN)

# ╔═╡ 47d445c5-d8d3-4c5f-8b17-4b6f302d86ce
fkz = fkz_7pts

# ╔═╡ 7be03d13-a754-43d1-9bc8-e1418f180c34
md"""
x coordinates in the air and in the slit
The center of the slit is at x=0
The slit is positioned in the air by Nm1 
"""

# ╔═╡ ffe03b6c-bbee-489e-91a0-507766775979
fxa(Nm1, Ns, dx, i) = (i - Nm1 - (Ns+1)/2)*dx

# ╔═╡ d4263ac0-5bd6-4d5a-ab66-99a35226f586
fxs(Ns, dx, i) = (i - (Ns+1)/2)*dx

# ╔═╡ a00d2be1-16cd-4f25-a6b4-9e87e89bfcf6
md"""
Solve the electric-flux equation
in the case of perfectly conducting resonant slit grating
using the 7-points t,x,z stencil.
The field is at normal incidence,
with a period equal to 1
    Ns is the number of sampling points in the slit

    The time period is 1
    The grating period is d, typically 0.5
    The slit width over d is Ns/Na, typically 1/6

    d, Nt, Na, Ns, dz = 0.5, 12*201, 6*201, 201, 1/(12*201)
    
"""

# ╔═╡ 2fd65b92-f674-46ed-a32e-ee1950107704
function compute_modes(p)
	d = p.d
	Nt = p.Nt
	Na = p.Na
	Ns = p.Ns
	dz = p.dz
    dt = 1/Nt    # distance between sampling points in t
    dx = d/Na
    ws = Ns*dx   # slit width
    
    
    Nm2 = (Na-Ns)÷2     # number of sampling points on the right part of the grating
    Nm1 = Na - Ns - Nm2 # number of sampling points on the left part of the grating
    
    xxa = fxa.(Nm1, Ns, dx, 1:Na) # x coordinates of the sampling points in the air
    xxs = fxs.(Ns, dx, 1:Ns)      # x coordinates of the sampling points in the slit
    
    vkxa = vcat(0:Na÷2,-(Na-1)÷2:1:-1)/d # kx wavenumbers in the air
    vkxs = (0:Ns-1)/(2*ws)               # kx wavenumbers in the slit
    
    vkza = fkz.(dt, dx, dz, vkxa) # kz wavenumbers in the air
    vkzs = fkz.(dt, dx, dz, vkxs) # kz wavenumbers in the slit

    """ control of the 7-points stencil """
    err = maximum(abs.((cos.(2pi*dt).-1)/dt^2/2 .- (cos.(2pi*vkxa*dx).-1)/dx^2/2 - (cos.(2pi*vkza*dz).-1)/dz^2/2))
    # println(err)
    @assert err*min(dt,dx,dz)^2 < 1e-13

    """
    matrix of order amplitudes in the air
    columns : orders
    rows : coordinates
    """
    Ma = Array{Complex{Float64}, 2}(undef, Na, Na)
    for i in 1:Na
        Ma[:,i] = exp.(1im*2*pi*vkxa[i]*xxa)
    end
    
    """
    matrix of mode amplitudes in the slit
    columns : modes
    rows : coordinates
    """
    Ms = Array{Float64, 2}(undef, Ns, Ns)
    for i in 1:2:Ns
        Ms[:,i] = cos.(2*pi*vkxs[i]*xxs)
    end
    for i in 2:2:Ns
        Ms[:,i] = sin.(2*pi*vkxs[i]*xxs)
    end
    
    Dm = Diagonal(vcat(ones(Nm1), zeros(Ns), ones(Nm2)))
    Ds = Diagonal(vcat(zeros(Nm1), ones(Ns), zeros(Nm2)))
    Rs = Ds[Nm1+1:Nm1+Ns,:]
    Rm1 = Dm[1:Nm1,:]
    Rm2 = Dm[Nm1+Ns+1:Na,:]
    
    DPhiap = Diagonal(exp.(1im*2pi*vkza*0.5*dz))
    DPhiam = Diagonal(exp.(-1im*2pi*vkza*0.5*dz))
    DPhisp = Diagonal(exp.(1im*2pi*vkzs*0.5*dz))
    DPhism = Diagonal(exp.(-1im*2pi*vkzs*0.5*dz))
    
    @assert maximum(abs.(DPhiap*DPhiam) .- 1) < 1e-14
    @assert maximum(abs.(DPhisp*DPhism) .- 1) < 1e-14
    
    Mc = Ms\(Rs*Ma)
    Mmr = DPhisp*Mc*DPhiap - DPhism*Mc*DPhiam
    Mmu = -DPhisp*Mc*DPhiam + DPhism*Mc*DPhiap
    Mad = Ma*(DPhiap - DPhiam) 
    Mad1 = Rm1*Mad
    Mad2 = Rm2*Mad
    Mr = vcat(Mad1, Mmr, Mad2)
    Mu = vcat(Mad1, Mmu, Mad2)
    Ua = vcat(1, zeros(Na-1))
    Ra = Mr\(Mu[:,1])
    Tsa = DPhisp*Mc*(DPhiam*Ua + DPhiap*Ra)
    Tsa_bis = DPhism*Mc*(DPhiap*Ua + DPhiam*Ra)
    @assert maximum(abs.(Tsa-Tsa_bis)) < 1e-14
    @assert abs(abs2(Ra[1]) + abs2(Tsa[1])*(Ns/Na) - 1) < 1e-14
    
    Mtm = vcat(Mad1, Rs*Ma*DPhiam, Mad2)
    Mtp = vcat(Mad1, Rs*Ma*DPhiap, Mad2)
    Mse = vcat(zeros(Nm1, Ns), Ms, zeros(Nm2, Ns))
    Mtms = Mtm\Mse
    Mtps = Mtp\Mse
    Rss = (Mtps*DPhism - Mtms*DPhisp)\((Mtms*DPhism - Mtps*DPhisp))
    Tas = Mtms*(DPhism + DPhisp*Rss)
    Tas_bis = Mtps*(DPhisp + DPhism*Rss)
    @assert maximum(abs.(Tas-Tas_bis)) < 1e-9
    
    @assert abs(abs2(Rss[1,1]) + abs2(Tas[1,1])*(Na/Ns) - 1) < 1e-9
    (Ra = Ra, Rss = Rss, Tas = Tas, Tsa = Tsa, Nt = Nt, Ns = Ns, Nm1 = Nm1, dt = dt, dx = dx, dz = dz, xxa = xxa, vkza = vkza)
end

# ╔═╡ 7d010001-2b4e-4aeb-9e5d-8ac81a161f7e
function compute_grid(p, limit=0)
	
	# Ra, Rss, Tas, Tsa, Nt, Ns, Nm1, dt, dx, dz, xxa, vkza
	
    raa = p.Ra[1]
    tas = p.Tas[1,1]
    tsa = p.Tsa[1]
    Us = vcat(1, zeros(p.Ns-1))

    As = -tas\raa
    Aa = p.Ra + p.Tas*(As*Us)
    if limit != 0
        for i in 2+limit:length(Aa)-limit
            Aa[i] = 0.0
        end
    end
    
    Nzz = p.Nt
    Ntt = p.Nt
    
    @assert isodd(p.Ns)
    
    i0 = p.Nm1 + (p.Ns+1)÷2
    
    @assert abs(p.xxa[i0]) < 1e-10
        
    zza = (0.5 .+ collect((Nzz-1):-1:0))*p.dz
    tta = collect(0:(Ntt-1))*p.dt
    Mzr = exp.(( 1im*2pi*zza) * transpose(p.vkza))*Aa
    Mzu = exp.((-1im*2pi*zza) * transpose(p.vkza)[1])
    Mz = Mzr + Mzu
    Mtz = real(Mz*transpose(exp.(-1im*2pi*tta)))
    (tta = tta, zza = zza, Nzz = Nzz, Mtz = Mtz, dt = p.dt, dz = p.dz)
end

# ╔═╡ 311c57eb-d00a-480b-80b8-f56d8c6ebbc2
# aev = 0.5457im
function analytic_grid(Ntt, Nzz, amplitude)
    Ntt = Ntt
    Nzz = Nzz
    dt = 1/Ntt
    dz = 1/Nzz
    
    zza = (0.5 .+ collect((Nzz-1):-1:0))*dz
    tta = collect(0:(Ntt-1))*dt
    Mzr = exp.(( 1im*2pi*zza) * im*sqrt(3.0))*2*amplitude
    Mzu = exp.((-1im*2pi*zza) * 1.0)
    Mz = Mzr + Mzu
    Mtz = real(Mz*transpose(exp.(-1im*2pi*tta)))
    (tta = tta, zza = zza, Nzz = Nzz, Mtz = Mtz, dt = dt, dz = dz)
end

# ╔═╡ 9797569a-7963-4d1d-87ce-90524f4d25bd
function itizfromtz(dt, dz, tmin, tmax, zmin, zmax)
    itmin = Int(1+round(tmin/dt))
    itmax = Int(round(tmax/dt))
    izmin = Int(1+round(zmin/dz))
    izmax = Int(round(zmax/dz))
    itmin, itmax, izmin, izmax
end

# ╔═╡ 6b5ffac5-96b0-4497-988e-67ce1bce4ec2
function zoom(p, tmin, tmax, zmin, zmax, levels)
    # en gros !
    itmin, itmax, izmin, izmax = itizfromtz(p.dt, p.dz, tmin, tmax, zmin, zmax)
    Plots.contour(p.tta[itmin:itmax], p.zza[p.Nzz:-1:1][izmin:izmax], p.Mtz[(p.Nzz:-1:1)[izmin:izmax],itmin:itmax],
                  levels=levels,
                  aspect_ratio=:equal,
                  show=true,
                  reuse=false)
end

# ╔═╡ 60aed1fc-7527-42cc-ba2b-1109179ee316
begin
    p01 = (d=0.5, Nt=12*201, Na=6*201, Ns=201, dz=1/(12*201))
    m01 = compute_modes(p01)
    g01 = compute_grid(m01, 1000)
end;

# ╔═╡ 847d6de0-f494-4001-9490-f180337531df
Plots.contour(g01.tta, g01.zza[g01.Nzz:-1:1], g01.Mtz[g01.Nzz:-1:1,:],
              levels= [
				  -0.81248,
				  -0.812,
				  -0.805,
				  -0.794,
				  -0.78,
				  -0.762,
				  -0.74,
				  -0.71,
				  -0.675,
				  -0.635],
               aspect_ratio=:equal)

# ╔═╡ fea6a6a6-9bec-4f10-ad72-f781287b8b58
zoom(g01, 0.3, 0.6, 0.01, 0.2, [
	-0.8149096,
	-0.8148,
	-0.812,
	-0.805,
	-0.794,
	-0.78,
	-0.762,
	-0.74,
	-0.71,
	-0.675,
	-0.635])

# ╔═╡ daa8b5bf-f08b-448a-91ee-89c8e554093a
zoom(g01, 0.415, 0.43, 0.115, 0.13, [
	-0.8149392,
	-0.8148,
	-0.8145,
	-0.8140,
	-0.8135])

# ╔═╡ cad320d8-fb15-4ee8-b4e1-4f9820251577
Plots.contour(g01.tta, g01.zza[g01.Nzz:-1:1], g01.Mtz[g01.Nzz:-1:1,:],
	levels=-0.826 .+ range(-0.5,0.5,step=0.1),
    aspect_ratio=:equal)

# ╔═╡ 4ae867e5-01ba-48ab-87f3-5912c1375bc1
zoom(g01, 0.3,0.6,0.01,0.2,-0.826 .+ range(-0.2,0.5,step=0.01))

# ╔═╡ 4c0d35aa-1622-48ef-ae70-a4e0a29f4785
ga = analytic_grid(2000,2000,0.5457im);

# ╔═╡ 849fbbf4-bcb6-40b2-bf03-9bcf3e8f5acc
Plots.contour(ga.tta, ga.zza[ga.Nzz:-1:1], ga.Mtz[ga.Nzz:-1:1,:],
    levels=-0.826 .+ range(-0.5,0.5,step=0.1),
    aspect_ratio=:equal)

# ╔═╡ d3e5fe26-50d3-4c5a-999c-159923a587e2
zoom(ga, 0.3,0.55,0.01,0.2,-0.826 .+ range(-0.2,0.5,step=0.01))

# ╔═╡ 1a584b16-a8c5-4d0d-aaf0-771e3a35bdd8
begin
    p02 = (d=0.5, Nt=12*401, Na=6*401, Ns=401, dz=1/(12*401))
    m02 = compute_modes(p02)
    g02 = compute_grid(m02)
end

# ╔═╡ afdcc8ad-e3a5-41d9-860a-d25ba20d0e8e
zoom(g02, 0.415, 0.43, 0.115, 0.13, [
	-0.8153994,
	-0.8150994,
	-0.8149497,
	-0.8148,
	-0.8145,
	-0.8140,
	-0.8135])

# ╔═╡ 4c89e9dc-ceb8-40d2-bedd-3874d829aa0e
zoom(g02, 0.415, 0.43, 0.115, 0.13, [-0.8149497])

# ╔═╡ a44f7fe9-79d5-47ac-8bff-4d6d8b52b1ef
# ╠═╡ disabled = true
#=╠═╡
let
    levels = [-0.8169, -0.8145,-0.8144, -0.812, -0.805, -0.794, -0.78, -0.762, -0.74, -0.71, -0.675, -0.635]
    
    Ns = 201
    Na = 6*Ns
    Nt = 2*Na
    dz = 1/Nt
    
    for i in 1:10
        # global dz, a
        dzo=dz
        a = compute_modes((d = 0.5, Nt = Nt, Na = Na, Ns = Na, dz = dz));
        Rss = a[2]
        rss = Rss[1,1]
        hs = (0-angle(rss)+pi)/2pi # 0.917
        Nzs = Int(round(hs/dz))
        dz = hs/Nzs
        println(dz-dzo)
    end
    
    Rss = a[2]
    rss = Rss[1,1]
    hs = (0-angle(rss)+pi)/2pi # 0.917
    Nzs = Int(round(hs/dz))
    dz = hs/Nzs

    g = compute_grid(a);

    Plots.contour(g.tta, g.zza[g.Nzz:-1:1], g.Mtz[g.Nzz:-1:1,:],
                  levels=levels,
                  aspect_ratio=:equal)

end
  ╠═╡ =#

# ╔═╡ 95527510-059b-48fd-ac85-769d3f106cc5
begin
    p03 = (d = 0.5, Nt = 14*201, Na = 6*201, Ns = 201, dz = 1/(18*201))
	m03 = compute_modes(p03)
	g03 = compute_grid(m03)
	l03 = [-128,-64,-32,-16,-8,-4,-1,0,1,4,8,16,32,64,128]*2e-6
end

# ╔═╡ 559d6453-4bdf-46fc-b841-ca5a6443e8b1
Plots.contour(g03.tta, g03.zza[g03.Nzz:-1:1], g03.Mtz[g03.Nzz:-1:1,:],
                  levels=l03,
                  aspect_ratio=:equal)

# ╔═╡ e525fae7-90e1-4d1d-880f-4ed143aed45a
zoom(g03, 0.415,0.43,0.115,0.13, -0.814943 .+ l03)

# ╔═╡ 6cf8584f-8937-45c3-ba04-b6316cceeef7
let
    itmin, itmax, izmin, izmax = itizfromtz(g03.dt, g03.dz, 0.420,0.426,0.115,0.13)
    Plots.plot(g03.Mtz[(g03.Nzz:-1:1)[izmin:izmax],(itmin+itmax)÷2])
end

# ╔═╡ 596d51f5-aaad-4010-a4c1-1a3837a6439e
Plots.plot(g03.Mtz[(g03.Nzz:-1:1)[435:455],1190])

# ╔═╡ Cell order:
# ╠═3db88667-117f-4485-b477-a6c886ce8c52
# ╠═a2228ad2-61da-4953-94a5-a5585b5c3a05
# ╠═8972d564-4a82-49b7-b9bc-b97812aa7f05
# ╟─7f4ee055-69aa-4603-a82c-d86e616721eb
# ╟─6a0001ef-699c-4669-8c0e-7cac66f40cc9
# ╠═fb6c7048-9566-4700-9c2b-9c1957225be6
# ╠═9cdbc786-3f96-4dce-bc9d-d3374ec7ab91
# ╠═47d445c5-d8d3-4c5f-8b17-4b6f302d86ce
# ╠═7be03d13-a754-43d1-9bc8-e1418f180c34
# ╠═ffe03b6c-bbee-489e-91a0-507766775979
# ╠═d4263ac0-5bd6-4d5a-ab66-99a35226f586
# ╟─a00d2be1-16cd-4f25-a6b4-9e87e89bfcf6
# ╠═2fd65b92-f674-46ed-a32e-ee1950107704
# ╠═7d010001-2b4e-4aeb-9e5d-8ac81a161f7e
# ╠═311c57eb-d00a-480b-80b8-f56d8c6ebbc2
# ╠═9797569a-7963-4d1d-87ce-90524f4d25bd
# ╠═6b5ffac5-96b0-4497-988e-67ce1bce4ec2
# ╠═60aed1fc-7527-42cc-ba2b-1109179ee316
# ╠═847d6de0-f494-4001-9490-f180337531df
# ╠═fea6a6a6-9bec-4f10-ad72-f781287b8b58
# ╠═daa8b5bf-f08b-448a-91ee-89c8e554093a
# ╠═cad320d8-fb15-4ee8-b4e1-4f9820251577
# ╠═4ae867e5-01ba-48ab-87f3-5912c1375bc1
# ╠═4c0d35aa-1622-48ef-ae70-a4e0a29f4785
# ╠═849fbbf4-bcb6-40b2-bf03-9bcf3e8f5acc
# ╠═d3e5fe26-50d3-4c5a-999c-159923a587e2
# ╠═1a584b16-a8c5-4d0d-aaf0-771e3a35bdd8
# ╠═afdcc8ad-e3a5-41d9-860a-d25ba20d0e8e
# ╠═4c89e9dc-ceb8-40d2-bedd-3874d829aa0e
# ╠═a44f7fe9-79d5-47ac-8bff-4d6d8b52b1ef
# ╠═95527510-059b-48fd-ac85-769d3f106cc5
# ╠═559d6453-4bdf-46fc-b841-ca5a6443e8b1
# ╠═e525fae7-90e1-4d1d-880f-4ed143aed45a
# ╠═6cf8584f-8937-45c3-ba04-b6316cceeef7
# ╠═596d51f5-aaad-4010-a4c1-1a3837a6439e
