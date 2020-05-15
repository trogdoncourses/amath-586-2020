
using LinearAlgebra, Plots, Printf, SparseArrays, LaTeXStrings

h = .1; m = convert(Int64,1/h)-1;
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
display(X)
display(Y)

vals = fill(0.,m,m)
U = map((x,y) -> x + m*(y-1), 1/h*X,1/h*Y)  # construct an enumeration of the grid for testing

tovec(U) = sparsevec(reshape(transpose(U[end:-1:1,:]), m^2, 1))[:]
tovec_a(U) = Array(reshape(transpose(U[end:-1:1,:]), m^2, 1))[:]
tomatrix(Uvec) = transpose(reshape(Uvec,m,m))[end:-1:1,:]

tovec(U) # The correct enumeration

U - tomatrix(tovec(U)) # test that tomatrix is the inverse of tovec

U

# not actually used -- the inverse of v_location
function bc_location(j::Int64,m1::Int64)
    row = convert(Int64,floor((j-1)/m1) + 1)
    col = j - (row-1)*m1
    (row,col)
end

function v_location(p,m) # p = (row, col)
    return (p[1]-1)*m+p[2] # basically the same function used to create the test grid U
end

function pos_gen(m)
    h0pos = map(t -> v_location(t,m), [(1,i) for i in 1:m]) # along the bottom boundary, h0
    h1pos = map(t -> v_location(t,m), [(m,i) for i in 1:m]) # along the top boundary, h1
    g0pos = map(t -> v_location(t,m), [(i,1) for i in 1:m]) # g0
    g1pos = map(t -> v_location(t,m), [(i,m) for i in 1:m]) # g1
    pos = [];
    append!(pos,h0pos)
    append!(pos,h1pos)
    append!(pos,g0pos)
    append!(pos,g1pos)
    pos
end    

pos = pos_gen(m)
pos[1:m] |> display # locations for h0(x,t), increasing x
pos[m+1:2m] |> display# locations for h1(x,t)
pos[2m+1:3m] |> display # locations for g0(y,t), increasing y
pos[3m+1:end] |> display # locations for g1(y,t)

function boundary(x::Vector{Float64},y::Vector{Float64},t::Float64,pos::Vector)
    bcvec = Array{Float64}([]);
    append!(bcvec,h0(x,t))
    append!(bcvec,h1(x,t))
    append!(bcvec,g0(y,t))
    append!(bcvec,g1(y,t))
    sparsevec(pos,bcvec)
    #Array(sparsevec(pos,bcvec))
end

h0 = (x,t) -> 1. .+ 0*x
h1 = (x,t) -> -1. .+ 0*x
g0 = (y,t) -> 2. .+ 0*y
g1 = (y,t) -> -2. .+ 0*y
U = boundary(x,y,0.,pos)
tomatrix(U) |> Array

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
A |> Array

bcfun = (x,y,t) -> t/(1 + t)*exp.(x.*y)*(sin(10*t)^2).+0.3
g0 = (y,t) -> bcfun(0.0,y,t)
g1 = (y,t) -> bcfun(1.0,y,t)
h0 = (y,t) -> bcfun(x,0.0,t)
h1 = (y,t) -> bcfun(x,1.0,t)
η = (x,y) -> 2*y.*x.*exp.(-30*(2*(x .- 0.5).^2 .+ (y .- 0.5).^2)).+0.3

function plot_heat(U,x,y,t,cl,width=800)
    l = @layout [a b]
    p1 = surface(x, y, U[end:-1:1,:], zaxis = [cl[1],cl[2]], clims= cl, aspectratio = .6, xlabel = L"x", ylabel = L"y", zlabel = L"u(x,y,t)")
    p2 = contour(x, y, U[end:-1:1,:], clims=cl, xaxis = [0,1], yaxis = [0,1], fill = true, aspectratio = 1, xlabel = L"x", ylabel = L"y")
    plot(p1, p2, layout = 2, size = (width, 7*width/10), title = @sprintf("t = %1.4f",t))
end

U⁰ = η(X,Y)
plot_heat(U⁰,x,y,0.0,(-0.1,1.5))

h = .01; m = convert(Int64,1/h)-1;
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
U⁰ = η(X,Y)
plot_heat(U⁰,x,y,0.0,(-0.1,1.5))

h = .01; k = h^2/4; T = 1;  #On the stability boundary
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
U⁰ = η(X,Y)
Ar = I + r*A;
pos = pos_gen(m)

U = tovec(U⁰)
t = 0
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U += r*boundary(x,y,t,pos)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(tomatrix(U),x,y,t,(-0.1,1.5)) |> IJulia.display
    end
end

h = .01; k = h^2/4*(1.001); T = .64;  #Just beyond the stability boundary
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
U⁰ = η(X,Y)
Ar = I + r*A;
pos = pos_gen(m)

U = tovec(U⁰)
t = 0
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = Ar*U
    U += r*boundary(x,y,t,pos)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(tomatrix(U),x,y,t,(-0.1,1.5)) |> IJulia.display
    end
end

function CG(A,b,x,eps,flag = false)
   r = b - A*x; p = r; n = 0
   while norm(r) > eps
        q = A*p
        a = (r'*r)/(p'*q)
        x = x + a*p
        r_old = r
        r = r - a*q
        b = (r'*r)/(r_old'*r_old)
        p = r + b*p 
        n += 1
    end
    if flag
        @printf("Iteration count = %i \n",n)
    end
    x
end

h = .01; k = .01; T = 2;
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(2*h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
Al = I - r*A; Ar = I + r*A;
pos = pos_gen(m);

U⁰ = η(X,Y)
plot_heat(U⁰,x,y,0.0,(-0.1,1.5))

U = tovec(U⁰)
t = 0
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    Uold = U
    U = Ar*U
    U += r*(boundary(x,y,t,pos)+boundary(x,y,t-k,pos))
    U = CG(Al,U,Uold,0.00001)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(tomatrix(U),x,y,t,(-0.1,1.5)) |> IJulia.display
    end
end

U⁰ = η(X,Y)
U = tovec(U⁰)
sol = CG(Al,Ar*U,U,0.00001,true)
Al*sol - Ar*U |> norm

off_diag = fill(1.0,m^2-1)
for i = m:m:m^2-1
    off_diag[i] = 0.
end
Ap = SymTridiagonal{Float64}(fill(-2,m^2),off_diag)
Ap[1:20,1:20] |> Array |> display
Ap = I - (r/4)*Ap;

## When A = Ap this will first "do" heat flow in the x-direction and then in y-direction.
function pre(A::SymTridiagonal{Float64,Array{Float64,1}},x::SparseVector{Float64,Int64})
    tovec(transpose(tomatrix(A\tovec_a(transpose(tomatrix(A\Array(x)))))))
end

function CG_pre(A::SparseMatrixCSC{Float64,Int64},M::Function,b::SparseVector{Float64,Int64},x::SparseVector{Float64,Int64},eps::Float64,flag = false)
   r = b - A*x; n = 0; z = M(r);  p = z;
   while norm(r) > eps
        q = A*p
        a = (r'*z)/(p'*q)
        x = x + a*p
        r_old = r
        z_old = z
        r = r - a*q
        z = M(r)
        b = (r'*z)/(r_old'*z_old)
        p = z + b*p 
        n += 1
    end
    if flag
        @printf("Iteration count = %i \n",n)
    end 
    x
end

U⁰ = η(X,Y)
U = tovec(U⁰)
@time sol = CG_pre(Al,x -> pre(Ap,x),Ar*U,U,0.0000001,true)
Al*sol - Ar*U |> norm

U⁰ = η(X,Y)
U = tovec(U⁰)
@time sol = CG(Al,Ar*U,U,0.0000001,true)
Al*sol - Ar*U |> norm

h = .01; k = .01; T = 1.;
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(2*h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
Al = I - r*A; Ar = I + r*A;
pos = pos_gen(m)

# Setup tridiagonal preconditioner
off_diag = fill(1.0,m^2-1)
for i = m:m:m^2-1
    off_diag[i] = 0.
end
Ap = SymTridiagonal{Float64}(fill(-2,m^2),off_diag)
Ap = I - (r/4)*Ap;

U = tovec(U⁰)
t = 0
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    Uold = U
    U = Ar*U
    U += r*(boundary(x,y,t,pos)+boundary(x,y,t-k,pos))
    U = CG_pre(Al, x-> pre(Ap,x),U,Uold,0.00001)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(tomatrix(U),x,y,t,(-0.1,1.5)) |> IJulia.display
    end
end

h = .01; k = .01; T = 1.;
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(2*h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
Al = I - r*A; Ar = I + r*A;
pos = pos_gen(m)

# Setup tridiagonal preconditioner
off_diag = fill(1.0,m^2-1)
for i = m:m:m^2-1
    off_diag[i] = 0.
end
Ap = SymTridiagonal{Float64}(fill(-2,m^2),off_diag)
Ap = I - (r/4)*Ap;
U⁰ = η(X,Y);

function with_prec()
    U = tovec(U⁰)
    t = 0
    for i = 2:n+1
        t += k
        Uold = U
        U = Ar*U
        U += r*(boundary(x,y,t,pos)+boundary(x,y,t-k,pos))
        U = CG_pre(Al, x-> pre(Ap,x),U,Uold,0.00001)
    end
end

function without_prec()
    U = tovec(U⁰)
    t = 0
    for i = 2:n+1
        t += k
        Uold = U
        U = Ar*U
        U += r*(boundary(x,y,t,pos)+boundary(x,y,t-k,pos))
        U = CG(Al,U,Uold,0.00001)
    end
end

@time with_prec()
@time without_prec()

h = .01; k = .001; T = 2.;
m = convert(Int64,1/h)-1;
n = convert(Int64,ceil(T/k))

A₀ = SymTridiagonal{Float64}(fill(-4.0,m),fill(1,m-1))
A = spzeros(m^2,m^2)
for j = 1:m
    A[m*(j-1)+1:m*j,m*(j-1)+1:m*j] += A₀
end
for j = 1:m-1
    A[m*(j-1)+1:m*j,m*(j)+1:m*(j+1)] += I
    A[m*(j)+1:m*(j+1),m*(j-1)+1:m*(j)] += I
end
r = k/(2*h^2)
x = Array(h:h:1-h)
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x))
Al = I - r*A; Ar = I + r*A;
pos = pos_gen(m)

# Setup tridiagonal preconditioner
off_diag = fill(1.0,m^2-1)
for i = m:m:m^2-1
    off_diag[i] = 0.
end
Ap = SymTridiagonal{Float64}(fill(-2,m^2),off_diag)
Ap = I - r*Ap;
U⁰ = η(X,Y);

plot()
anim = Animation()
U = tovec(U⁰)
t = 0
plot_heat(tomatrix(U),x,y,0.0,(-0.1,1.5))
frame(anim)
fr = 1000 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    Uold = U
    U = Ar*U
    U += r*(boundary(x,y,t,pos)+boundary(x,y,t-k,pos))
    U = CG_pre(Al, x-> pre(Ap,x),U,Uold,0.00001)
    if mod(i-1,tb) ≈ 0.0
        #IJulia.clear_output(true)
        plot_heat(tomatrix(U),x,y,t,(-0.1,1.5))
        frame(anim)
    end
end
gif(anim,"heat_2D_CN.gif")
