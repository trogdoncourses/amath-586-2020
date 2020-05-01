
using LinearAlgebra, Plots, Printf

h = 0.01;
m = convert(Int64,1/h)-1;
k = 0.01;
T = 10;

A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
r = k/(2*h^2)
Al = I - r*A
Ar = I + r*A;

g0 = t -> t.^2/(1 .+ t.^2)*sin.(4*t)
g1 = t -> 0.
η = x -> exp.(-20*(x .-1/2).^2)
println(g0(0.))
println(η(0.))

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1-h
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U[1] += r*(g0(t)+g0(t-k))
    U[2] += r*(g1(t)+g1(t-k))
    U = Al\U
    if mod(i-1,tb) ≈ 0.0 # working with i*k makes it so our framerate doesn't change with k
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter)
        frame(anim)
    end
end
gif(anim,"heat_CN.gif")

h = 0.01;
m = convert(Int64,1/h)-1;
k = 0.01;
T = 10;
A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
r = k/(2*h^2)
Al = I - r*A
Ar = I + r*A;

n = convert(Int64,ceil(T/k))
x = h:h:1-h
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U[1] += r*(g0(t)+g0(t-k))
    U[2] += r*(g1(t)+g1(t-k))
    U = Al\U
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

function CN_heat(η,g0,g1,T,k,h)
    m = convert(Int64,1/h)-1;
    A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
    r = k/(2*h^2)
    Al = I - r*A
    Ar = I + r*A;
    n = convert(Int64,ceil(T/k))
    x = h:h:1-h
    U = η(x)
    t = 0.0

    for i = 2:n+1
        t += k
        U = Ar*U
        U[1] += r*(g0(t)+g0(t-k))
        U[2] += r*(g1(t)+g1(t-k))
        U = Al\U
    end
    U
end

@time CN_heat(η,g0,g1,10.,0.0001,0.0001);

h = 0.01;
m = convert(Int64,1/h)-1;
k = h^2/2;  # boundary of stability region
T = 10;
A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
r = k/(h^2) # note that this is doubled
# Al = I - r*A # no need for this
Ar = I + r*A;

n = convert(Int64,ceil(T/k))
x = h:h:1-h
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U[1] += r*(g0(t))
    U[2] += r*(g1(t))
    # U = Al\U remove this step
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

h = 0.01;
m = convert(Int64,1/h)-1;
k = h^2*0.50002;  # just outside of stability region
T = 100;
A = SymTridiagonal(fill(-2.0,m),fill(1.0,m-1))
r = k/(h^2) # note that this is doubled
# Al = I - r*A # no need for this
Ar = I + r*A;

n = convert(Int64,ceil(T/k))
x = h:h:1-h
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display

fr = 10 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U[1] += r*(g0(t))
    U[2] += r*(g1(t))
    # U = Al\U remove this step
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),g1(t)], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end
