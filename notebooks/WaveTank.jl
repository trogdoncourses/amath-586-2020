
using LinearAlgebra, Plots, Printf, SparseArrays, LaTeXStrings

a = 0.; b = 10
h = 0.01
m = convert(Int64,b/h) #increase by 1 over Dirichlet conditions
T = 4
κ = .1;
k = h

B = Tridiagonal(fill(1.0,m-1),fill(-2.0,m),fill(1.0,m-1))
B[end,end-1] = 2
B *= 1/(h^2)
B1 = I - κ*(k/2)*B
B2 = I + κ*(k/2)*B;
r = κ*k/(2*h^2)

g0 = t -> 2*sin(2*pi*t).^2
η = x -> 0.

n = convert(Int64,ceil(T/k))
x = h:h:b |> Array
U = map(η,x)
t = 0.0
anim = Animation()
plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    # diffuse with Crank-Nicolson
    U = (B2*U)
    U[1] += r*(g0(t)+g0(t-k))
    U = B1\U
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
        plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"no_advection.gif")

a = 0.; b = 10
h = 0.01/2
m = convert(Int64,b/h) #increase by 1 over Dirichlet conditions
T = 4
α = 3.
κ = .1;
k = h/(2α)

B = Tridiagonal(fill(1.0,m-1),fill(-2.0,m),fill(1.0,m-1))
B[end,end-1] = 2
B *= 1/(h^2)
B1 = I - κ*(k/2)*B
B2 = I + κ*(k/2)*B;
r = κ*k/(2*h^2)

A = Tridiagonal(fill(-1.0,m-1),fill(0.0,m),fill(1.0,m-1))
A[end,end-1] = 0
A *= -α/(2h);
A1 = I + k*A + (k^2*α^2)/(2)*B; # Lax-Wendroff matrix

g0 = t -> 2*sin(2*pi*t).^2
η = x -> 0.

n = convert(Int64,ceil(T/k))
x = h:h:b |> Array
U = map(η,x)
t = 0.0
anim = Animation()
plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    
    ## O(k) splitting method
    # diffuse with Crank-Nicolson
    U = (B2*U)
    U[1] += r*(g0(t)+g0(t-k))
    U = B1\U
    #advect with Lax-Wendroff
    U = A1*U
    U[1] += (α*k/(2h) + (α^2*k^2/(2h^2)))*g0(t-k)
    
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
        plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"advection_diffusion.gif")

a = 0.; b = 10
h = 0.01
m = convert(Int64,b/h) #increase by 1 over Dirichlet conditions
T = 6
α = 3.
κ = 0.01;
k = h/(2α)

B = Tridiagonal(fill(1.0,m-1),fill(-2.0,m),fill(1.0,m-1))
B[end,end-1] = 2
B *= 1/(h^2)
B1 = I - κ*(k/2)*B
B2 = I + κ*(k/2)*B;
r = κ*k/(2*h^2)

A = Tridiagonal(fill(-1.0,m-1),fill(0.0,m),fill(1.0,m-1))
A[end,end-1] = 0
A *= -α/(2h);
A1 = I + k*A + α^2*k^2/(2)*B; # Lax-Wendroff matrix

g0 = t -> t < .5 ? 2*sin(2*pi*t).^2 : 0.0
η = x -> 0.

n = convert(Int64,ceil(T/k))
x = h:h:b |> Array
U = map(η,x)
t = 0.0

anim = Animation()
plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    
    ## O(k) splitting method
    # diffuse with Crank-Nicolson
    U = (B2*U)
    U[1] += r*(g0(t)+g0(t-k))
    U = B1\U
    #advect with Lax-Wendroff
    U = A1*U
    U[1] += (α*k/(2h) + (α^2*k^2/(2h^2)))*g0(t-k)
    
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
        plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"linear_wavetank.gif")

a = 0.; b = 5
h = 0.01
m = convert(Int64,b/h) #increase by 1 over Dirichlet conditions
T = 6
α = 2.
κ = 0.01;
k = h/4

B = Tridiagonal(fill(1.0,m-1),fill(-2.0,m),fill(1.0,m-1))
B[end,end-1] = 2
B *= 1/(h^2)
B1 = I - κ*(k/2)*B
B2 = I + κ*(k/2)*B;
r = κ*k/(2*h^2)

A = Tridiagonal(fill(-1.0,m-1),fill(0.0,m),fill(1.0,m-1))
A[end,end-1] = 0
A *= -α/(2h);

g0 = t -> t < 1. ? 2*sin(2*pi*t)^2 : 0.0
η = x -> 0.

n = convert(Int64,ceil(T/k))
x = h:h:b |> Array
U = map(η,x)
t = 0.0
anim = Animation()
plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    # diffusion
    U = (B2*U)
    U[1] += r*(g0(t)+g0(t-k))
    U = B1\U
    # nonlinear advection
    AU = A*U
    U1 = U[1]
    U2 = U[2]
    U = U + k*U.*AU + k^2*U.*(AU.^2) + (α^2*k^2/2)*(U.^2).*(B*U)
    g = g0(t-k)
    H = (α*k/(2h))*U1*g - (α^2*k^2/(2h^2))*U1*U2*g + (α^2*k^2/(4h^2))*U1*g^2 + (α^2*k^2/(2h^2))*U1^2*g
    U[1] += H
    
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [a,b], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
        plot!([a,b],[g0(t),U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"nonlinear_wavetank.gif")

ay = 0.; by = 1
hy = 0.01
my = convert(Int64,(by-ay)/hy) + 1 #increase by 2 over Dirichlet conditions
T = 2
α = 2.
κ = 0.01;
k = hy/4

By = Tridiagonal(fill(1.0,my-1),fill(-2.0,my),fill(1.0,my-1))
By[end,end-1] = 2
By[1,2] = 2
By *= 1/(hy^2)
B1y = I - κ*(k/2)*By
B2y = I + κ*(k/2)*By;
ry = κ*k/(2*hy^2)

Ay = Tridiagonal(fill(-1.0,my-1),fill(0.0,my),fill(1.0,my-1))
Ay[1,2] = 0
Ay[end,end-1] = 0
Ay *= -α/(2hy);

η = x -> exp(-25(x-.5)^2)

n = convert(Int64,ceil(T/k))
y = ay:hy:by |> Array
U = map(η,y)
t = 0.0
anim = Animation()
plot(y, U, xaxis = [ay,by], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
plot!([ay,by],[U[1],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = B1y\(B2y*U)

    AU = Ay*U
    U = U + k*U.*AU + k^2*U.*(AU.^2) + (α^2*k^2/2)*(U.^2).*(By*U)
    
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(y, U, xaxis = [ay,by], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t), fill = (-2,:lightblue))
        plot!([ay,by],[U[1],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"burgers_neumann.gif")

ax = 0.0; bx = 5.0;
ay = 0.0; by = 1.0;
hx = 0.1/10
hy = 0.1/2
x = ax+hx:hx:bx |> Array
y = ay:hy:by |> Array
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x));
Y = Y[end:-1:1,:]
αx = 2.0
αy = .1
κ = 0.01
k = min(hx,hy)/(4*max(αx,αy))

my = convert(Int64,(by-ay)/hy) + 1 #increase by 2 over Dirichlet conditions
By = Tridiagonal(fill(1.0,my-1),fill(-2.0,my),fill(1.0,my-1))
By[end,end-1] = 2
By[1,2] = 2
By *= 1/(hy^2)
B1y = I - κ*(k/2)*By
B2y = I + κ*(k/2)*By;
ry = κ*k/(2*hy^2)

Ay = Tridiagonal(fill(-1.0,my-1),fill(0.0,my),fill(1.0,my-1))
Ay[1,2] = 0
Ay[end,end-1] = 0
Ay *= -αy/(2hy);

mx = convert(Int64,(bx-ax)/hx) #increase by 1 over Dirichlet conditions
Bx = Tridiagonal(fill(1.0,mx-1),fill(-2.0,mx),fill(1.0,mx-1))
Bx[end,end-1] = 2
Bx *= 1/(hx^2)
B1x = I - κ*(k/2)*Bx
B2x = I + κ*(k/2)*Bx;
rx = κ*k/(2*hx^2)

Ax = Tridiagonal(fill(-1.0,mx-1),fill(0.0,mx),fill(1.0,mx-1))
Ax[end,end-1] = 0
Ax *= -αx/(2hx);

function plot_wave(U,x,y,t,cl,width=800)
    l = @layout [a;b]
    p1 = surface(x, y, U, camera = (10,70), zaxis = [cl[1],cl[2]], clims= cl, aspectratio = 1.5, xlabel = L"x", ylabel = L"y", zlabel = L"u(x,y,t)")
    p2 = contour(x, y, U, clims=cl, fill = true, aspectratio = 1, xlabel = L"x", ylabel = L"y")
    plot(p1, p2, layout = (2,1), size = (width, 7*width/10), title = @sprintf("t = %1.4f",t))
end

η = (x,y) -> 0*exp(-25*(x-2)^2-50*(y-.5)^2)
g0 = (y,t) -> 2*sin(2*pi*t).^2 .*exp(-25*(y-.5)^2)

U = map(η,X,Y)
T = 10
n = convert(Int64,ceil(T/k))
t = 0.0

cl = (-1,2)
plot_wave(U,x,y,t,cl) |> IJulia.display


anim = Animation()
frame(anim)
fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    # diffuse in y
    U = B1y\(B2y*U)
    # advect in y
    AU = Ay*U
    U = U + k*U.*AU + k^2*U.*(AU.^2) + (αy^2*k^2/2)*(U.^2).*(By*U)
   
    U = U' |> Array 
    # diffuse in x
    U = (B2x*U)
    U[1,:] += map(yy -> rx*(g0(yy,t)+g0(yy,t-k)),y)
    U = B1x\U
    # advect in x    
    AU = Ax*U
    U1 = U[1,:]
    U2 = U[2,:]
    U = U + k*U.*AU + k^2*U.*(AU.^2) + (αx^2*k^2/2)*(U.^2).*(Bx*U)
    g = map(yy -> g0(yy,t-k), y)
    H = (αx*k/(2h))*U1.*g - (αx^2*k^2/(2h^2))*U1.*U2.*g + (αx^2*k^2/(4h^2))*U1.*g.^2 + (αx^2*k^2/(2h^2))*U1.^2 .*g
    U[1,:] += H  
    
    U = U' |> Array
    
    
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_wave(U,x,y,t,cl) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"viscous_Burgers_2D.gif")
