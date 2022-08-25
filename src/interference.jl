# see also 2019-whirls/src/slit-modal.jl

import Pkg; Pkg.activate("env")
# Pkg.add("Plots")
# Pkg.add("Contour")
# Pkg.add("NLsolve")
# Pkg.add("DifferentialEquations")
# include("src/interference.jl")

using LinearAlgebra
import Plots
import Contour
using NLsolve
using DifferentialEquations

function itizfromtz(dt, dz, tmin, tmax, zmin, zmax)
    itmin = Int(1+round(tmin/dt))
    itmax = Int(round(tmax/dt))
    izmin = Int(1+round(zmin/dz))
    izmax = Int(round(zmax/dz))
    itmin, itmax, izmin, izmax
end


function zoom(dt, dz, tmin, tmax, zmin, zmax, p)
    # en gros !
    itmin, itmax, izmin, izmax = itizfromtz(dt, dz, tmin, tmax, zmin, zmax)
    Plots.contour(tta[itmin:itmax], zza[Nzz:-1:1][izmin:izmax], Mtz[(Nzz:-1:1)[izmin:izmax],itmin:itmax],
                  levels=p,
                  aspect_ratio=:equal,
                  show=true,
                  reuse=false)
end

# amplitude = 0.5457im
function analytic(Ntt, Nzz, imaga)
    Ntt = Ntt
    Nzz = Nzz
    dt = 1/Ntt
    dz = 1/Nzz
    
    zza = (0.5 .+ collect((Nzz-1):-1:0))*dz
    tta = collect(0:(Ntt-1))*dt
    Mtz = Array{Float64,2}(undef, Nzz, Ntt)
    for it in 1:Ntt
        t = tta[it]
        for iz in 1:Nzz
            z = zza[iz]
            Mtz[iz, it] = phie(t, z, imaga)
        end
    end
    tta, zza, Nzz, Mtz, dt, dz
end


tta, zza, Nzz, Mtz, dt, dz = analytic(2000,2000,IMAGA)
Plots.closeall()
Plots.contour(tta, zza[Nzz:-1:1], Mtz[Nzz:-1:1,:],
              levels=-0.826 .+ range(-0.5,0.5,step=0.1),
              aspect_ratio=:equal)

zoom(dt, dz, 0.3,0.55,0.01,0.2,-0.826 .+ range(-0.2,0.5,step=0.01))

IMAGA = 0.5457
phie(t, z, a) = cospi(2*(t+z)) +
                2*a*exp(-2pi*sqrt(3)*z)*sinpi(2*t)

hy(t, z, a) = 2pi*(-sinpi(2*(t+z)) +
              2*a*exp(-2pi*sqrt(3)*z)*cospi(2*t))
dx(t, z, a) = 2pi*(sinpi(2*(t+z)) +
              2*sqrt(3)*a*exp(-2pi*sqrt(3)*z)*sinpi(2*t))
phie(t, z) = phie(t, z, IMAGA)
hy(t, z) = hy(t, z, IMAGA)
dx(t, z) = dx(t, z, IMAGA)
(phie(t+dt,z)-phie(t,z))/dt, hy(t, z)
-(phie(t,z+dz)-phie(t,z))/dz, dx(t, z)




# recherche du centre de tourbillon Hy=0 Dx=0

function f!(F, x)
    F[1] = hy(x[1],x[2])
    F[2] = dx(x[1],x[2])
end

"""
function j!(J, x)
    J[1, 1] = x[2]^3-7
    J[1, 2] = 3*x[2]^2*(x[1]+3)
    u = exp(x[1])*cos(x[2]*exp(x[1])-1)
    J[2, 1] = x[2]*u
    J[2, 2] = u
end
"""

t = 5/12 # Obtenu à la main
z = 0.112313448736877236

nz1(z) = -log(-sin(2pi*(z+5/12))/(sqrt(3)*IMAGA))/(2pi*sqrt(3))

nz2(z) = 0.5-asin(-sqrt(3)*IMAGA*exp(-2pi*sqrt(3)*z))/2pi-5/12
function zval(z)
    for i in 1:100
        z = nz2(z)
    end
    z
end
z = zval(0.123)

# aussi nlsolve(f!,[5/12, 0.12313448736877236])

function spaghetti!(du,u,p,t)
 du[1] = hy(u[1],u[2])
 du[2] = dx(u[1],u[2])
end


function spag(t0, z0, span)
    u0 = [t0 ; z0]
    tspan = (0.0, span)
    prob = ODEProblem(spaghetti!, u0, tspan)
    sol = solve(prob, Rosenbrock23(), reltol=1e-10, abstol=1e-10)

    npts = length(sol.t)
    sss = sol.t
    ttt = Vector{Float64}(undef, npts)
    zzz = Vector{Float64}(undef, npts)
    for i in 1:npts
        ttt[i] = sol.u[i][1]
        zzz[i] = sol.u[i][2]
    end
    sss, ttt, zzz
end

sss,ttt,zzz = spag(0.4, 0.25, 1.5) 

zoom(dt, dz, 0.2,0.6,0.01,0.25,-0.826 .+ range(-0.2,0.5,step=0.01))
Plots.plot!(ttt,zzz)
