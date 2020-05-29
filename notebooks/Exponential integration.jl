
using FFTW, Plots, LinearAlgebra, LaTeXStrings, SparseArrays, Printf
include("./amath586.jl")
using .amath586
import Base: diff

L = 11
m = 10
X = -L .+ 2*L*(0:m-1)/m

rnd = xx -> map(x -> round(x, digits = 5),xx)
mfftshift = x -> circshift(fftshift(x), isodd(length(x)) ? 1 : 0)
mfft = x -> fftshift(fft(fftshift(x)))
mifft = x -> mfftshift(ifft(mfftshift(x)))
mgrid = (n,L) -> -L .+ 2*L*(0:n-1)/n

f = x -> exp.(0*1im*pi*x/L)
abs.(mfft(f(mgrid(m,L)))) |> rnd

f = x -> exp.(-1*1im*pi*x/L)
abs.(mfft(f(mgrid(m,L)))) |> rnd

f = x -> exp.(1*1im*pi*x/L)
abs.(mfft(f(mgrid(m,L)))) |> rnd

u = 1.0:10 |> Array
u - mifft(mfft(u)) |> rnd

m = 6
L = 10
f = x -> exp.(1*1im*pi*x/L)
out = mfft(f(mgrid(m,L))) |> rnd

m = 5
L = 10
f = x -> exp.(1*1im*pi*x/L)
out = mfft(f(mgrid(m,L))) |> rnd

struct trig_interp
    L::Float64
    c::Vector{Complex{Float64}}
end

function (tr::trig_interp)(x)
    m = length(tr.c)
    mm = convert(Int64,floor( m/2 ))
    σ = isodd(m) ? 1im*pi/m : 0. # if n is odd we need to rotate coefs
    ex = exp.(-1im*pi*mm*x/tr.L + mm*σ)
    ex1 = exp.(1im*pi*x/tr.L-σ)
    sum = tr.c[1]*ex
    for i = 2:length(tr.c)
        ex  =  ex.*ex1
        sum += tr.c[i]*ex
    end
    return sum/m
end    

L = 10.
n = 51
X = mgrid(n,L)
f = x -> exp.(-cos.(pi*x/L))
fm = trig_interp(L,mfft(f(X)));

f(.1)-fm(.1)

L = 10
f = x -> exp.(-10*cos.(pi*x/L))

data = []
ms = []
for m = 11:5:150
    X = mgrid(m,L)
    fn = trig_interp(L,mfft(f(X)))
    append!(data,abs(fn(.1)-f(.1)))
    append!(ms,m)
end

plot(ms,data,yaxis = :log, label = "", grid = true, framestyle = :box)
plot!(ms,data,yaxis = :log, seriestype = :scatter, label = "Errors at x = 0.1")

diffvec = (L,m,j) -> ((-floor(m/2):1:floor((m-1)/2))*(1im*pi/L)).^j

function diff(tr::trig_interp,j=1)
    return trig_interp(tr.L,diffvec(L,length(tr.c),j).*tr.c)
end

L = 10.
m = 51
X = mgrid(m,L)
f = x -> exp.(-cos.(pi*x/L))
df = x -> f(x).*(pi/L)*sin.(pi*x/L)
fm = trig_interp(L,mfft(f(X)));
dfm = diff(fm);

dfm(.1) - df(.1)

anim = Animation()
η = x -> exp.(-cos.(4*pi*x/L))

L = 50.
m = 2^10
γ = 1.
X = mgrid(m,L)
c = mfft(η(X))
D3 = γ*diffvec(L,m,3)
U = mifft(c)
cl = [-2,5]
k = 0.1
t = 0.0


plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)

for i = 2:200
    t += k
    U = mifft(exp.(-D3*t).*c)
    IJulia.clear_output(true)        
    plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
    frame(anim)
end
gif(anim,"linear_kdv_expcos.gif")

anim = Animation()
η = x -> 2*exp.(-3x.^2)

L = 50.
m = 2^10
γ = 2.
X = mgrid(m,L)
c = mfft(η(X))
D3 = γ*diffvec(L,m,3)
U = mifft(c)
cl = [-2,5]
k = 0.01
t = 0.0


plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)

for i = 2:200
    t += k
    U = mifft(exp.(-D3*t).*c)
    IJulia.clear_output(true)        
    plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
    frame(anim)
end
gif(anim,"linear_kdv_gauss.gif")

function rk4(F,k,t,c)
    f1 = k*F(c,t)
    f2 = k*F(c + .5*f1, t + .5*k)
    f3 = k*F(c + .5*f2, t + .5*k)
    f4 = k*F(c + f3, t + k)
    return c + 1/6.0*(f1 + 2.0*f2 + 2.0*f3 + f4)
end

anim = Animation()
T = 4
η = x -> exp.(-cos.(4*pi*x/L))

L = 50.
m = 2^10
k = 0.0001

γ = 1.
X = mgrid(m,L)
c = mfft(η(X))
U = mifft(c)
cl = [-1,6]

D = diffvec(L,m,1)
Am = -γ*diffvec(L,m,3)

F = (v,τ) -> -3*exp.(-Am*τ).*(D.*mfft(mifft(exp.(Am*τ).*v).^2))

n = convert(Int64,ceil(T/k))
t = 0.0
plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    c = exp.(Am*k).*rk4(F,k,0.0,c)
    t += k
    U = mifft(c)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)        
        plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
        frame(anim)
    end
end
gif(anim,"kdv_expcos.gif")


anim = Animation()
T = 4
η = x -> 2*exp.(-3x.^2)

L = 50.
m = 2^10
k = 0.0001

γ = 1.
X = mgrid(m,L)
c = mfft(η(X))
U = mifft(c)
cl = [-1,6]

D = diffvec(L,m,1)
Am = -γ*diffvec(L,m,3)

F = (v,τ) -> -3*exp.(-Am*τ).*(D.*mfft(mifft(exp.(Am*τ).*v).^2))

n = convert(Int64,ceil(T/k))
t = 0.0
plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    c = exp.(Am*k).*rk4(F,k,0.0,c)
    t += k
    U = mifft(c)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)        
        plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
        frame(anim)
    end
end
gif(anim,"kdv_gauss.gif")

anim = Animation()
T = 4
A = 3.0
η = x -> A*sech.(sqrt(A/2)*x).^2

L = 50.
m = 2^10
k = 0.0001

γ = 1.
X = mgrid(m,L)
c = mfft(η(X))
U = mifft(c)
cl = [-4,4]

D = diffvec(L,m,1)
Am = -γ*diffvec(L,m,3)

F = (v,τ) -> -3*exp.(-Am*τ).*(D.*mfft(mifft(exp.(Am*τ).*v).^2))

n = convert(Int64,ceil(T/k))
t = 0.0
plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    c = exp.(Am*k).*rk4(F,k,0.0,c)
    t += k
    U = mifft(c)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)        
        plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
        frame(anim)
    end
end
gif(anim,"kdv_sech.gif")


anim = Animation()
T = 4
A = 3.0
η = x -> -A*sech.(sqrt(A/2)*x).^2

L = 50.
m = 2^10
k = 0.0001

γ = 1.
X = mgrid(m,L)
c = mfft(η(X))
U = mifft(c)
cl = [-4,4]

D = diffvec(L,m,1)
Am = -γ*diffvec(L,m,3)

F = (v,τ) -> -3*exp.(-Am*τ).*(D.*mfft(mifft(exp.(Am*τ).*v).^2))

n = convert(Int64,ceil(T/k))
t = 0.0
plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display
frame(anim)
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    c = exp.(Am*k).*rk4(F,k,0.0,c)
    t += k
    U = mifft(c)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)        
        plot(X, U |> real, xaxis = [-L,L], yaxis = cl, lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (cl[1],:lightblue)) |> IJulia.display  
        frame(anim)
    end
end
gif(anim,"kdv_neg_sech.gif")

