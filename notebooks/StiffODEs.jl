
using Plots, LaTeXStrings, LinearAlgebra
include("./amath586.jl")
using .amath586

## Small κ
κ = .1;  
f = u -> [u[1] - u[1]*u[2] + (κ/2)*u[1]*u[2]*u[3], u[1]*u[2] - u[2] + (κ/2)*u[1]*u[2]*u[3], -κ*u[1]*u[2]*u[3] ]

# Forward Euler
k = 0.001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t,U[1,:],label=L"u_1(t)")
plot!(t,U[2,:],label=L"u_2(t)")
plot!(t,U[3,:],label=L"u_3(t)",yaxis=[0,7],grid=true)

# Forward Euler
k = 0.001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,5.] # 5 instead of 3
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t,U[1,:],label=L"u_1(t)")
plot!(t,U[2,:],label=L"u_2(t)")
plot!(t,U[3,:],label=L"u_3(t)",yaxis=[0,7],grid=true)

κ = 200.;  # Increase κ
f = u -> [u[1] - u[1]*u[2] + (κ/2)*u[1]*u[2]*u[3], u[1]*u[2] - u[2] + (κ/2)*u[1]*u[2]*u[3], -κ*u[1]*u[2]*u[3] ]

# Forward Euler
k = 0.001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[0,7],grid=true)

# Forward Euler with smaller step size
k = 0.0001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,3.]
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot(t[1:10:end],U[1,1:10:end],label=L"u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"u_3(t)",yaxis=[0,7],grid=true)

# Forward Euler
k = 0.001
n = convert(Int64,50/k)
U = zeros(3,n+1)
t = zeros(n+1)
U[:,1] = [5.,.8,0.]
t[1] = 0.
for i = 2:n+1
    U[:,i] = U[:,i-1] + k*f(U[:,i-1])
    t[i] = t[i-1] + k
end

println("Maximum value of u_1(t):  ",maximum(U[1,:]))
plot!(t[1:10:end],U[1,1:10:end],label=L"\tilde u_1(t)")
plot!(t[1:10:end],U[2,1:10:end],label=L"\tilde u_2(t)")
plot!(t[1:10:end],U[3,1:10:end],label=L"\tilde u_3(t)",yaxis=[0,7],grid=true)
