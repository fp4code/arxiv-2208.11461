using Pkg
Pkg.activate("../env")

"""
Pkg.add("TaylorSeries")
Pkg.add("ForwardDiff")
Pkg.add("NLsolve")
Pkg.add("Plots")
Pkg.add("Contour")
"""

using TaylorSeries
#using NLSolve
#using ForwardDiff

vt, vz = set_variables("t z", order=2)

phie(a, t, z) = sin(z+t) - a*exp(-sqrt(3)*z)*cos(t)


fz(a, z) = cos(z)/sqrt(3) - sin(z) + a*exp(-sqrt(3)*z)
dfz(a, z) = -sin(z)/sqrt(3) - cos(z) - a*sqrt(3)*exp(-sqrt(3)*z)
newz(a, z) = z - fz(a, z)/dfz(a, z)

a = 2*0.54561

function goodz(n, a, z)
    for i in 1:n
        z = newz(a, z)
    end
    z
end

t0 = 2pi/6
z0 = goodz(20, a, 0.12312*2pi)

phier(a, t, z) = phie(a, t0+t, z0+z)

att, atz, azz = phier(a, vt, vz).coeffs[3].coeffs

(-sin(z0+t0) + a*exp(-sqrt(3)*z0)*cos(t0))/2 - att
(-sin(z0+t0) - 3*a*exp(-sqrt(3)*z0)*cos(t0))/2 - azz
(-sin(z0+t0) - sqrt(3)*a*exp(-sqrt(3)*z0)*sin(t0)) - atz
(-sin(z0+t0) - 3*a*exp(-sqrt(3)*z0)*cos(t0)) - atz

cos(t0)/2 - 1/4



S = sin(z0+t0)/2
S = (sqrt(3)*cos(z0) + sin(z0))/4
E = a*exp(-sqrt(3)*z0)/4
-S + E - att
-S - 3*E - azz
-S - 3*E - atz/2


import Plots
import Contour
using LinearAlgebra
import Plots

Ntt = 101
Nzz = 101

tt = range(-pi,pi,length=Ntt)
zz = range(0,2pi,length=Nzz)

m = Array{Float64,2}(undef, Nzz, Ntt)

for it in 1:Ntt
    for iz in 1:Nzz
        m[iz, it] = phie(a, tt[it], zz[iz])
    end
end

ttr = range(t0-pi,t0+pi,length=Ntt)
zzr = range(-z0,2pi-z0,length=Nzz)

m = Array{Float64,2}(undef, Nzz, Ntt)

for it in 1:Ntt
    for iz in 1:Nzz
        m[iz, it] = phie(a, t0+ttr[it], z0+zzr[iz])
    end
end

Plots.contour(ttr/2pi, zzr/2pi, m)

ttr = range(-z0,z0,length=Ntt)
zzr = range(-z0,z0,length=Nzz)

m = Array{Float64,2}(undef, Nzz, Ntt)
mr = Array{Float64,2}(undef, Nzz, Ntt)

for it in 1:Ntt
    for iz in 1:Nzz
        m[iz, it] = phie(a, t0+ttr[it], z0+zzr[iz])
        mr[iz, it] = att*ttr[it]^2 + atz*ttr[it]*zzr[iz] + azz*zzr[iz]^2
    end
end

Plots.contour(ttr/2pi, zzr/2pi, m)
Plots.contour(ttr/2pi, zzr/2pi, mr,
              levels=range(-2,2,step=0.05),
              aspect_ratio=:equal)

zoft(pm,t) = (-1 + pm*sqrt(1 - att/azz))*t

Plots.plot!(ttr/2pi, zoft(1,ttr)/2pi)
Plots.plot!(ttr/2pi, zoft(-1,ttr)/2pi)

(S-E)/(S+3E) - att/azz

Czoftp = (-1 + sqrt(1 - att/azz))
Czoftm = (-1 - sqrt(1 - att/azz))
Ctofzp = (-1 - sqrt(1 - att/azz))*azz/att
Ctofzm = (-1 + sqrt(1 - att/azz))*azz/att

Czoftp*Ctofzp - 1.0
Czoftm*Ctofzm - 1.0


aa = att/azz
A = [aa 1; -1 -1]*2*azz

# eta -> 0 -inf

# https://mathworld.wolfram.com/MatrixExponential.html
delta = sqrt(-(aa + 1)^2+4)

V = [(1-aa - sqrt((1 - aa)^2 - 4+0im))/2 (1-aa + sqrt((1 - aa)^2 - 4+0im))/2; 1 1]
LL = [(-1 + aa + sqrt((1 - aa)^2 - 4+0im))/2 0 ;0 (-1 + aa - sqrt((1 - aa)^2 - 4+0im))/2]*2*azz


expA(S, E, eta) = exp([(S-E)/(S+3E) 1; -1 -1]*2*(-S-3E)*eta)

expA(S, E, -2.0)


function expAe(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    ddi = sqrt(aa^2 +2aa - 3 + 0im)
    h = 2*(-S-3E)*eta
    AA = ((ddi + aa + 1)/2)*exp((ddi + aa - 1)*h/2) - ((-ddi + aa + 1)/2) * exp((-ddi + aa - 1)*h/2)
    AB = exp((ddi + aa - 1)*h/2) - exp((-ddi + aa - 1)*h/2)
    AC = ((ddi - aa - 1)/2)*exp((ddi + aa - 1)*h/2) - ((-ddi - aa - 1)/2) * exp((-ddi + aa - 1)*h/2)
    [AA AB;-AB AC]/ddi
end


function expA01(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    ddi = dd*1im
    h = 2*(-S-3E)*eta
    AA = ((ddi + aa + 1)/2)*exp(ddi*h/2) - ((-ddi + aa + 1)/2) * exp(-ddi*h/2)
    AB = exp(ddi*h/2) - exp(-ddi*h/2)
    AC = ((ddi - aa - 1)/2)*exp(ddi*h/2) - ((-ddi - aa - 1)/2) * exp(-ddi*h/2)
    [AA AB;-AB AC]*exp((aa-1)*h/2)/ddi
end

function expA02(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    ddi = dd*1im
    h = 2*(-S-3E)*eta
    AA = ((ddi + aa + 1)/2im)*exp(ddi*h/2) - ((-ddi + aa + 1)/2im) * exp(-ddi*h/2)
    AB = (exp(ddi*h/2) - exp(-ddi*h/2))/1im
    AC = ((ddi - aa - 1)/2im)*exp(ddi*h/2) - ((-ddi - aa - 1)/2im) * exp(-ddi*h/2)
    [AA AB;-AB AC]*exp((aa-1)*h/2)/dd
end

function expA03(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    ddi = dd*1im
    h = 2*(-S-3E)*eta
    AA = ((ddi + aa + 1)/2im)*exp(ddi*h/2) - ((-ddi + aa + 1)/2im) * exp(-ddi*h/2)
    AB = 2*sin(dd*h/2)
    AC = ((ddi - aa - 1)/2im)*exp(ddi*h/2) - ((-ddi - aa - 1)/2im) * exp(-ddi*h/2)
    [AA AB;-AB AC]*exp((aa-1)*h/2)/dd
end

function expA04(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    ddi = dd*1im
    h = 2*(-S-3E)*eta
    AA = ((ddi + aa + 1)/2im)*exp(ddi*h/2) - ((-ddi + aa + 1)/2im) * exp(-ddi*h/2)
    AA = dd*cos(dd*h/2) + (aa+1)*sin(dd*h/2)
    AB = 2*sin(dd*h/2)
    AC = ((ddi - aa - 1)/2im)*exp(ddi*h/2) - ((-ddi - aa - 1)/2im) * exp(-ddi*h/2)
    [AA AB;-AB AC]*exp((aa-1)*h/2)/dd
end

function expA05(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    h = 2*(-S-3E)*eta
    AA = dd*cos(dd*h/2) + (aa + 1)*sin(dd*h/2)
    AB = 2*sin(dd*h/2)
    AC = dd*cos(dd*h/2) - (aa + 1)*sin(dd*h/2)
    [AA AB;-AB AC]*exp((aa-1)*h/2)/dd
end

function expA06(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    h = 2*(-S-3E)*eta
    AA = dd*cos(dd*h/2) + (aa + 1)*sin(dd*h/2)
    AB = 2*sin(dd*h/2)
    AC = dd*cos(dd*h/2) - (aa + 1)*sin(dd*h/2)
    ([1 0; 0 1]*dd*cos(dd*h/2) + [aa+1 2; -2 -(aa+1)]*sin(dd*h/2))*exp((aa-1)*h/2)/dd
end

function expA07(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    h = (-S-3E)*eta
    AA = dd*cos(dd*h) + (aa + 1)*sin(dd*h)
    AB = 2*sin(dd*h)
    AC = dd*cos(dd*h) - (aa + 1)*sin(dd*h)
    ([1 0; 0 1]*dd*cos(dd*h) + [aa+1 2; -2 -(aa+1)]*sin(dd*h))*exp((aa-1)*h)/dd
end

function expA08(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    AA = dd*cos(dd*azz*eta) + (aa + 1)*sin(dd*azz*eta)
    AB = 2*sin(dd*azz*eta)
    AC = dd*cos(dd*azz*eta) - (aa + 1)*sin(dd*azz*eta)
    ([1 0; 0 1]*dd*cos(dd*azz*eta) + [aa+1 2; -2 -(aa+1)]*sin(dd*azz*eta))*exp((aa-1)*azz*eta)/dd
end

function expA09(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2aa + 3)
    azzdd = azz*dd
    ([1 0; 0 1]*cos(azzdd*eta) + [(aa+1)/2 1; -1 -(aa+1)/2]*(2/dd)*sin(azzdd*eta))*exp((aa-1)*azz*eta)
end

function expA10(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2*aa + 3)
    azzdd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    ([1 0; 0 1]*cos(azzdd*eta) + [(aa+1)/2 1; -1 -(aa+1)/2]*(2/dd)*sin(-azzdd*eta))*exp(-(aa/dd-1/dd)*azzdd*eta)
end

function expA11(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2*aa + 3)
    azzdd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    ([1 0; 0 1]*cos(azzdd*eta) + [(aa+1)/2 1; -1 -(aa+1)/2]*(2/dd)*sin(-azzdd*eta))*exp((1-aa)*(1/dd)*azzdd*eta)
end

function expA12(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2*aa + 3)
    azzdd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    ([1 0; 0 1]*cos(azzdd*eta) + [(aa+1)/2 1; -1 -(aa+1)/2]*(2/dd)*sin(-azzdd*eta))*exp(-(1-aa)*azz*eta)
end

function expA13(S, E, eta)
    att = -S + E
    azz = -S - 3E
    aa = att/azz
    dd = sqrt(-aa^2 - 2*aa + 3)
    azzdd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    ([1 0; 0 1]*cos(azzdd*eta) + [(aa+1)/2 1; -1 -(aa+1)/2]*(2/dd)*sin(-azzdd*eta))*exp((att - azz)*eta)
end

function expA14(S, E, eta)
    att = -S + E
    azz = -S - 3E
    dd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    ([1 0; 0 1]*cos(dd*eta) + [(att+azz)/2 azz; -azz -(att+azz)/2]*(2/dd)*sin(dd*eta))*exp((att - azz)*eta)
end

function expA15(S, E, eta)
    att = -S + E
    azz = -S - 3E
    dd = sqrt(32*E^2 + 16*S*E)
    ([1 0; 0 1]*cos(dd*eta) + [(att+azz)/2 azz; -azz -(att+azz)/2]*(2/dd)*sin(dd*eta))*exp((att - azz)*eta)
end

function expA16(S, E, eta)
    att = -S + E
    azz = -S - 3E
    dd = sqrt(-att^2 - 2*att*azz + 3*azz^2)
    dd = 4*sqrt(E*(S+2*E))
    ([1 0; 0 1]*cos(dd*eta) + [(att+azz)/2 azz; -azz -(att+azz)/2]*(2/dd)*sin(dd*eta))*exp((att - azz)*eta)
end


function expA17(S, E, eta)
    dd = 4*sqrt(E*(S+2E))
    ([1 0; 0 1]*cos(dd*eta) + [-(S+E) -(S+3E); S+3E S+E]*(2/dd)*sin(dd*eta))*exp(4*E*eta)
end

function expAaz0_01(a, z0, eta)
    e = a*exp(-sqrt(3)*z0)
    s = sqrt(3)*cos(z0) + sin(z0)
    dd = sqrt(e*(s+2e))
    ([1 0; 0 1]*cos(dd*eta) + [-(s+e) -(s+3e); s+3e s+e]*(1/(2*dd))*sin(dd*eta))*exp(e*eta)
end

expAaz0_01(a, z0, -2.0) - expA(S, E, -2.0)

function expAaz0_02(z0, eta)
    e = -(1/sqrt(3))*cos(z0) + sin(z0)
    s = sqrt(3)*cos(z0) + sin(z0)
    dd = sqrt(e*(s+2e))
    ([1 0; 0 1]*cos(dd*eta) + [-(s+e) -(s+3e); s+3e s+e]*(1/(2*dd))*sin(dd*eta))*exp(e*eta)
end

expAaz0_02(z0, -2.0) - expA(S, E, -2.0)

function expAaz0_03(z0, eta)
    ad = sin(z0) - (1/sqrt(3))*cos(z0)
    s1 = sin(z0) + (1/sqrt(3))*cos(z0)
    s3 = 2*sin(z0)
    dd = sqrt((sin(z0) - (1/sqrt(3))*cos(z0))*(3*sin(z0) + (1/sqrt(3))*cos(z0)))
    ([1 0; 0 1]*cos(dd*eta) + [-s1 -s3; s3 s1]*(1/dd)*sin(dd*eta))*exp(ad*eta)
end

expAaz0_03(z0, -2.0) - expA(S, E, -2.0)

function expAaz0_04(z0, eta)
    ad = sin(z0) - (1/sqrt(3))*cos(z0)
    dd = sqrt(-sqrt(3)*sin(2*z0) - 5*cos(2*z0) + 4)/sqrt(3)
    s1 = (sin(z0) + (1/sqrt(3))*cos(z0))/dd
    s3 = 2*sin(z0)/dd
    ([1 0; 0 1]*cos(dd*eta) + [-s1 -s3; s3 s1]*sin(dd*eta))*exp(ad*eta)
end

expAaz0_04(z0, -2.0) - expA(S, E, -2.0)



function tzmat(z0)
    ad = sin(z0) - (1/sqrt(3))*cos(z0)
    dd = sqrt(-sqrt(3)*sin(2*z0) - 5*cos(2*z0) + 4)/sqrt(3)
    s1 = (sin(z0) + (1/sqrt(3))*cos(z0))/dd
    s3 = 2*sin(z0)/dd
    [1 0; 0 1], [-s1 -s3; s3 s1], dd, ad
end

function ploplo(mimi, mama)

    Plots.contour(ttr/2pi, zzr/2pi, mr,
                  levels=range(-2,2,step=0.05),
                  aspect_ratio=:equal)
        
    Plots.plot!(ttr/2pi, zoft(1,ttr)/2pi)
    Plots.plot!(ttr/2pi, zoft(-1,ttr)/2pi)
    
    
    mc, ms, dd, ad = tzmat(z0)

    eta = range(mimi,mama, step=-0.1)
    
    
    for xx in range(0, 0.1, length=11)
        ts = -xx*2pi
        zs = xx*2pi
        ad = 0
        vc = [ts ; zs]
        vs = ms*vc
        tt = exp.(ad.*eta).*(cos.(dd.*eta)*vc[1] + sin.(dd.*eta)*vs[1])
        zz = exp.(ad.*eta).*(cos.(dd.*eta)*vc[2] + sin.(dd.*eta)*vs[2])
        
        Plots.plot!(tt/2pi, zz/2pi)
    end
end

"""
t = u + v
z = u - v

u = (cos.(dd.*eta) u_s + sin.(dd.*eta)*(s3-s1)*v_s)exp.(ad.*eta)
v = (cos.(dd.*eta) v_s + sin.(dd.*eta)*(-s3-s1)*u_s)exp.(ad.*eta)

Or (s3-s1)*(s3+s1) = 1

u' = k u
v' = (1/k) v

k = sqrt(s3 + s1) = `1/sqrt(s3-s1)

t = t' (k + 1/k)/2 + z' (-k + 1/k)/2
z = z' (-k + 1/k)/2 + z' (k + 1/k)/2

v = v/z pour z' = 0
"""

k = sqrt(s3 + s1)      #  1.7214597300856138
(k + 1/k)/2            #  1.1511810392768584
(-k + 1/k)/2           # -0.5702786908087554
v = (-k+1/k)/(k + 1/k) # -0.49538575719331635

mrp = Array{Float64,2}(undef, Nzz, Ntt)
for it in 1:Ntt
    for iz in 1:Nzz
        t = ttr[it] * (k + 1/k)/2 + zzr[iz] * (-k + 1/k)/2
        z = ttr[it] * (-k + 1/k)/2 + zzr[iz] * (k + 1/k)/2
        mrp[iz, it] = att*t^2 + atz*t*z + azz*z^2
    end
end

Plots.contour(ttr/2pi, zzr/2pi, mrp,
              levels=range(-2,2,step=0.05),
              aspect_ratio=:equal)
mc, ms, dd, ad = tzmat(z0)

mimi = -5*(2pi/dd)
mama = 2*(2pi/dd)
eta = range(mimi,mama, step=0.1)
    
    
for xx in range(0, 2pi*ad/dd, length=11)
    ts = -exp(xx)*2pi/exp(2pi*ad/dd)*0.1
    zs = exp(xx)*2pi/exp(2pi*ad/dd)*0.1
    vc = [ts ; zs]
    vs = ms*vc
    tt = exp.(ad.*eta).*(cos.(dd.*eta)*vc[1] + sin.(dd.*eta)*vs[1])
    zz = exp.(ad.*eta).*(cos.(dd.*eta)*vc[2] + sin.(dd.*eta)*vs[2])
    ttp = tt * (k + 1/k)/2 + zz * (k - 1/k)/2
    zzp = tt * (k - 1/k)/2 + zz * (k + 1/k)/2
   
    Plots.plot!(ttp/2pi, zzp/2pi, legend = false)
end

Plots.plot!(0,0, xlims = (-0.1, 0.1), ylims = (-0.1, 0.1)) 
