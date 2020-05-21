
using LinearAlgebra, Plots, Printf, SparseArrays
include("./amath586.jl")
using .amath586

h = 0.01
m = convert(Int64,1/h)-1;
k = 0.1
T = 10.
a = 1.0;

# Initial condition
η = x -> exp.(-20*(x .-1/2).^2)
# Initial condition chosen so that u(x,t) = sin(2*pi*(x - t)), if a = 1
# η = x -> sin.(2*pi*x)
# u = (x,t) -> sin.(2*pi*(x.-t))

A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)); #increase dim by 1
A = sparse(A) # Need to convert A₀ to a new data type to allow new entries
A[1,end] = -1
A[end,1]  = 1
A *= -a/(2h)
A |> Array;

k = 2*h^2
B = I + k*A;

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B*U
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_periodic_FE.gif")

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
t += k
Uold = U
U = B*U

for i = 3:n+1
    t += k
    Utemp = U
    U = Uold + 2*k*A*U
    Uold = Utemp
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_periodic_leap.gif")

λ = eigvals(k*A |> Array)
plot(real(λ),imag(λ), seriestype = :scatter, xaxis = [-1,1])

a*k/h # bound on the modulus of the eigenvalues

# Leapfrog
α = [1,0,-1]
β = [0,2,0]
convergence_stability(α,β)

h = 0.01
m = convert(Int64,1/h)-1;
k = h/a
T = 10.
a = 1.0;
η = x -> exp.(-20*(x .-1/2).^2)
A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)); #increase dim by 1
A = sparse(A) # Need to convert A₀ to a new data type to allow new entries
A[1,end] = -1
A[end,1]  = 1
A *= -a/(2h);

a*k/h # bound on the modulus of the eigenvalues

λ = eigvals(k*A |> Array)
plot(real(λ),imag(λ), seriestype = :scatter, xaxis = [-1,1])

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
t += k
Uold = U
U = B*U

for i = 3:n+1
    t += k
    Utemp = U
    U = Uold + 2*k*A*U
    Uold = Utemp
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_periodic_leap_stable.gif")

h = 0.01
m = convert(Int64,1/h)-1;
k = h/a+0.00001
T = 2.
a = 1.0;
η = x -> exp.(-20*(x .-1/2).^2)
A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)); #increase dim by 1
A = sparse(A) # Need to convert A₀ to a new data type to allow new entries
A[1,end] = -1
A[end,1]  = 1
A *= -a/(2h);

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
t += k
Uold = U
U = B*U

for i = 3:n+1
    t += k
    Utemp = U
    U = Uold + 2*k*A*U
    Uold = Utemp
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[end],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_periodic_leap_unstable.gif")

# Initial and boundary conditions
η = x -> exp.(-20*(x .-1/2).^2)
g0(t) = sin(4*t)

h = 0.01
m = convert(Int64,1/h)-1;
k = h^2/4
T = 10.
a = 1.0;

A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m));
A = sparse(A)
A[end,end-2:end] = [1,-4,3];
A *= -a/(2h);

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
# If solution is known
# plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    
    U = U + k*(A*U)
    U[1] += a*k/(2h)*g0(t-k)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        # If solution is known
        # plot!(x, u(x,t), xaxis = [0,1], yaxis = [-1,2],lw=1,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_dirichlet_FE.gif")


h = 0.01
m = convert(Int64,1/h)-1;
k = h^2/4
T = 5.
a = 1.0;

A = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)); #increase dim by 1
A = sparse(A)
A[end,end-2:end] = [1,-4,3];
A *= -a/(2h);

anim = Animation()
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter)
frame(anim)

fr = 50 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
t += k
Uold = U
U = U + k*(A*U)
U[1] += a*k/(2h)*g0(t-k)

for i = 3:n+1
    t += k
    Utemp = U
    U = Uold + 2*k*A*U
    U[1] += a*k/(h)*g0(t-k)
    Uold = Utemp
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
IJulia.clear_output(true)
gif(anim,"advection_dirichlet_leap.gif")


λ = eigvals(k*A |> Array)
plot(real(λ),imag(λ), seriestype = :scatter, xaxis = [-1,1])

λ = eigvals(k*A |> Array)
plot(real(λ),imag(λ), seriestype = :scatter, xaxis = [-.001,.001])
