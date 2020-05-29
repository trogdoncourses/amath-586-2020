
using LinearAlgebra, Plots, Printf

h = 0.01;
m = convert(Int64,1/h)-1;
k = 0.01;
T = 10;

A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
r = 1im*k/(2*h^2)
Al = I - r*A
Ar = I + r*A;

g0 = t -> t.^2/(1 .+ t.^2)*sin.(4*t)
g1 = t -> 0.
η = x -> exp.(-20*(x .-1/2).^2)

plot()
anim = Animation()
n = convert(Int64,ceil(T/abs(k)))
x = h:h:1-h
U = η(x)
t = 0.0
plot(x, U |> real, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("Re u(x,t), t = %1.2f",t))
plot!(x, U |> imag, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("Im u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),g1(t)] |> real, label = "Real BCs", seriestype = :scatter)
plot!([0,1],[g0(t),g1(t)] |> imag, label = "Imag BCs", seriestype = :scatter)


frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
     t += k
     U = Ar*U
     U[1] += r*(g0(t)+g0(t-k))
     U[2] += r*(g1(t)+g1(t-k))
     U = Al\U
     if mod(i-1,tb) ≈ 0.0
         plot(x, U |> real, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("Re u(x,t), t = %1.2f",t))
         plot!(x, U |> imag, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("Im u(x,t), t = %1.2f",t))
         plot!([0,1],[g0(t),g1(t)] |> real, label = "Real BCs", seriestype = :scatter)
         plot!([0,1],[g0(t),g1(t)] |> imag, label = "Imag BCs", seriestype = :scatter)
         frame(anim)
     end
 end
gif(anim,"LS_CN.gif")

function TR_LS(η,g0,g1,T,k,h)
    m = convert(Int64,1/h)-1;
    A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
    r = 1im*k/(2*h^2)
    Al = I - r*A
    Ar = I + r*A;
    n = convert(Int64,ceil(T/k))
    x = h:h:1-h
    U = η(x)
    t = 0.0

    @inbounds for i = 2:n+1
        t += k
        U = Ar*U
        U[1] += r*(g0(t)+g0(t-k))
        U[2] += r*(g1(t)+g1(t-k))
        U = Al\U
    end
    U
end

@time TR_LS(η,g0,g1,.25,.25/10^3,0.001); # 10^3 time steps, h = 1/10^3
