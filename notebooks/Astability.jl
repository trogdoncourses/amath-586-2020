
using Plots, LaTeXStrings, LinearAlgebra
include("./amath586.jl")
using .amath586

κ = 200.;  # Increase κ
f = u -> [u[1] - u[1]*u[2] + (κ/2)*u[1]*u[2]*u[3], u[1]*u[2] - u[2] + (κ/2)*u[1]*u[2]*u[3], -κ*u[1]*u[2]*u[3] ]
Df = u -> [1.0-u[2]+(κ/2)*u[2]*u[3] -u[1]+(κ/2)*u[1]*u[3] (κ/2)*u[1]*u[2];
    u[2]+(κ/2)*u[2]*u[3] u[1]-1.0+(κ/2)*u[1]*u[3] (κ/2)*u[1]*u[2];
    -κ*u[2]*u[3] -κ*u[1]*u[3] -κ*u[1]*u[2]  ]

g = (u,Uⁿ) -> u - Uⁿ - (k/2)*(f(u)+f(Uⁿ))
Dg = u -> I - (k/2)*Df(u)

k = 0.001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    # u -> g(u,U[:,i-1]) produces the function we want the root of
    U[:,i] = run_newton(u->g(u,U[:,i-1]),Dg,U[:,i-1],k^2/10,10) 
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[0,7],grid=true)

plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[-2,2],grid=true)

eigvals(Df([5.,.8,3.]))

g = (u,Uⁿ) -> u - Uⁿ - (k)*(f(u))
Dg = u -> I - (k)*Df(u)

k = 0.002
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    # u->g(u,U[:,i-1]) produces the function we want the root of
    U[:,i] = run_newton(u->g(u,U[:,i-1]),Dg,U[:,i-1],k^2/10,20) 
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[0,7],grid=true)

plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[-2,2],grid=true)

# BDF r = 3
α = [11,-18,9,-2]
β = [6,0,0,0]
convergence_stability(α,β)

c = .9; a = c -> (1-2c)/(2*(1-c)); b = c -> 1/(2*(1-c));
R = z -> (1 + b(c)*z/(1-c*z))/(1-b(c)*z)
xrange = [-3,3]; yrange = [-3,3]
contourf(xrange[1]:0.01(1+rand()/10):xrange[2],yrange[1]:0.01(1+rand()/10):yrange[2],(x,y)-> sign(1-(abs(R(x+1im*y)))),colorbar=false)

g1 = (u,Uⁿ,Uᵐ,Uˡ) -> 11*u - 18*Uⁿ + 9*Uᵐ - 2*Uˡ - (6k)*f(u)  # for BDF r = 3
g2 = (u,Uⁿ) -> u - Uⁿ - c*k*f(u) # first stage of DIRK
g3 = (u,U,Uⁿ) -> u - Uⁿ - k*(b(c)*f(U)+a(c)*f(u)) #second stage of DIRK
Dg1 = u -> 11*I-(6k)*Df(u)
Dg2 = u -> I - (k/3)*Df(u)
Dg3 = u -> I - (k/4)*Df(u)

k = 0.0001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
max_iter = 30
for i = 2:3
    t[i] = t[i-1] + k
    Ustar = run_newton(u -> g2(u,U[:,i-1]),Dg2,U[:,i-1],1e-14,max_iter) # first RK stage
    U[:,i] = run_newton(u -> g3(u,Ustar,U[:,i-1]),Dg3,Ustar,1e-14,max_iter) # second RK stage
end
for i = 4:n+1
    t[i] = t[i-1] + k
    U[:,i] = run_newton(u -> g1(u,U[:,i-1],U[:,i-2],U[:,i-3]),Dg1,U[:,i-1],1e-14,max_iter) # use BDF r = 3
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[0,7],grid=true)

Utrue = [3.663299609166387,0.02453535459551013,1.04e-322]; # this was the output

k = 0.0004
N = 7
approx = zeros(N)
for l = 1:N
    k = k/2
    n = convert(Int64,ceil(10/k))
    U = zeros(3,n+1)
    t = zeros(n+1)
    U[:,1] = [5.,.8,3.]
    t[1] = 0.
    max_iter = 10
    for i = 2:3
        t[i] = t[i-1] + k
        Ustar = run_newton(u -> g2(u,U[:,i-1]),Dg2,U[:,i-1],1e-14,100)
        U[:,i] = run_newton(u -> g3(u,Ustar,U[:,i-1]),Dg3,Ustar,1e-14,100)
    end
    for i = 4:n+1
        t[i] = t[i-1] + k
        U[:,i] = run_newton(u -> g1(u,U[:,i-1],U[:,i-2],U[:,i-3]),Dg1,Ustar,1e-14,100)
    end
    approx[l] = maximum(abs.(U[:,end]-Utrue))
    if l > 1
        println("Error reduction ratio  ",approx[l-1]/approx[l])
    end
end

R = (r,w₀,z) -> ChebyshevT(r,w₀ + z*ChebyshevT(r,w₀)/DChebyshevT(r,w₀))/ChebyshevT(r,w₀)

w₀ = 1; r = 3;
xrange = [-80,.1]; yrange = [-10,10]
contourf(xrange[1]:0.1(1+rand()/10):xrange[2],yrange[1]:0.1(1+rand()/10):yrange[2],(x,y)-> sign(1-(abs(R(r,w₀,x+1im*y)))),colorbar=false)

r = 6; w₀ = 1 + 0.05/r^2
xrange = [-80,.1]; yrange = [-10,10]
contourf(xrange[1]:0.1(1+rand()/10):xrange[2],yrange[1]:0.1(1+rand()/10):yrange[2],(x,y)-> sign(1-(abs(R(r,w₀,x+1im*y)))),colorbar=false)
