
using LinearAlgebra, Plots, Printf, SparseArrays

h = 0.01
m = convert(Int64,1/h)-1;
k = 0.1
T = 10.
A₀ = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)) #increase dim by 1
a = 1.0;

A₀ |> Array;

# Initial condition
η = x -> exp.(-20*(x .-1/2).^2)
# Initial condition chosen so that u(x,t) = sin(2*pi*(x - t)), if a = 1
# η = x -> sin.(2*pi*x)
# u = (x,t) -> sin.(2*pi*(x.-t))

A = sparse(A₀) # Need to convert A₀ to a new data type to allow new entries
A[1,end] = -1
A[end,1]  = 1
A *= -a/(2h)
A |> Array;

k = h^2
B = I + k*A;

n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
# If solution is known
# plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        # If solution is known
        # plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    if mod(i-1,tb) ≈ 0.0
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
        frame(anim)
    end
end
gif(anim,"advection_periodic.gif")

# Initial condition
η = x -> exp.(-20*(x .-1/2).^2)
# Initial condition chosen so that u(x,t) = sin(2*pi*(x - t)), if a = 1
#η = x -> sin.(2*pi*x)
#u = (x,t) -> sin.(2*pi*(x.-t))

h = 0.001
m = convert(Int64,1/h)-1;
T = 10.
A₀ = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)) #increase dim by 1
a = 1.0;

A = sparse(A₀) # Need to convert A₀ to a new data type to allow new entries
A[1,end] = -1
A[end,1]  = 1
A *= -a/(2h)
C₀ = sparse(Tridiagonal(fill(0.5,m),fill(0.0,m+1),fill(0.5,m)))
C₀[1,end] = .5
C₀[end,1] = .5
k = h # need this to avoid too much numerical dissipation
B = sparse(C₀ + k*A);

n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
# If solution is known
# plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        # If solution is known
        # plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    if mod(i-1,tb) ≈ 0.0
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
        frame(anim)
    end
end
gif(anim,"advection_periodic_LF.gif")

# Initial and boundary conditions
η = x -> exp.(-20*(x .-1/2).^2)
g0(t) = sin(4*t)+1.

h = 0.01
a = 1.0;
m = convert(Int64,1/h)-1;
k = h/(a)
T = 10.
A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)); #increase dim by 1
A = sparse(A)
A[end,end-2:end] = [1,-4,3];

C = Tridiagonal(fill(0.5,m),fill(0.0,m+1),fill(0.5,m))
C = sparse(C)
C[end,:] *= 0.0
C[end,end] = 1.0
B = sparse(C - (a*k/(2h))*A);

n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    U[1] += (1/2 + a*k/(2h))*g0(t-k)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    U[1] += (1 + a*k/h)*g0(t)/2
    if mod(i-1,tb) ≈ 0.0
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)
        frame(anim)
    end
end
gif(anim,"advection_LF.gif")

# Initial and boundary conditions
η = x -> exp.(-40*(x .-1/2).^2)
g0(t) = sin(4*t)^2.
# Initial condition chosen so that u(x,t) = sin(2*pi*(x - t)), if a = 1
# η = x -> sin.(2*pi*x)
# g0(t) = sin.(-2*pi*t)
# u = (x,t) -> sin.(2*pi*(x.-t))

h = 0.01
a = 1.0;
m = convert(Int64,1/h)-1;
k = h/(a)
T = 10.
A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)) |> sparse; 
D = Tridiagonal(fill(1.0,m),fill(-2.0,m+1),fill(1.0,m)) |> sparse;
vec = [1,-4,3]
A[end,end-2:end] = vec; # same as Lax-Friedrichs
# Construct a backward second-order approximation of the second derivative
D[end,:] *= 0.0
D[end,end-3:end-1] += (vec[1]/2)*[-.5,0,.5]
D[end,end-2:end] += (vec[2]/2)*[-.5,0,.5]
D[end,end-2:end] += vec[3]*vec/4;

n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
# If solution is known
# plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)

fr = 1000 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = U - (a*k/(2h))*(A*U) + (a^2*k^2/(2h^2))*(D*U)
    U[1] += (a*k/(2h) + (a^2*k^2/(2h^2)))*g0(t-k)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        # If solution is known
        # plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
    end
end

# Initial and boundary conditions chosen so that u(x,t) is the solution

F = x -> (exp.(-(x .- 1).^2) .+ 1).*sin.(2*pi*x)
η = x -> F(x)
g0(t) = F(-t)
u = (x,t) -> F(x.-t)

h = 0.2
p = 7
out = fill(0.0,p)
for i = 1:p
    h = h/2
    a = 1.0;
    m = convert(Int64,1/h)-1;
    k = h/(a)
    T = 10.
    A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)) |> sparse; 
    D = Tridiagonal(fill(1.0,m),fill(-2.0,m+1),fill(1.0,m)) |> sparse;
    vec = [1,-4,3]
    A[end,end-2:end] = vec; # same as Lax-Friedrichs
    # Construct a backward second-order approximation of the second derivative
    D[end,:] *= 0.0
    D[end,end-3:end-1] += (vec[1]/2)*[-.5,0,.5]
    D[end,end-2:end] += (vec[2]/2)*[-.5,0,.5]
    D[end,end-2:end] += vec[3]*vec/4;
    
    n = convert(Int64,ceil(T/k))
    x = h:h:1 # include right end point
    U = η(x)
    t = 0.0

    for i = 2:n+1
        t += k
        U += -(a*k/(2h))*(A*U) + (a^2*k^2/(2h^2))*(D*U)
        U[1] += (a*k/(2h) + (a^2*k^2/(2h^2)))*g0(t-k)
    end
    out[i] = maximum(abs.(U - u(x,T)))
end
out[1:end-1]./out[2:end]
