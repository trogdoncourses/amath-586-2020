
using LinearAlgebra, Plots, Printf, LaTeXStrings

h = 0.01
κ = 0.3
m = convert(Int64,1/h)+1; #increase by 2 over Dirichlet conditions
k = 0.001;
T = 1;
B = Tridiagonal(fill(1.0,m-1),fill(-2.0,m),fill(1.0,m-1))
B[1,2] = 2
B[end,end-1] = 2;
B *= κ/(h^2);

B1 = I - (k/4)*B
B2 = I + (k/4)*B
B3 = I - (k/3)*B;
function TRBDF2_heat(U)
    Ustar = B1\(B2*U)
    return B3\((4/3)*Ustar-(1/3)*U)
end

η = x -> exp.(-20*(x .-1/2).^2)

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))
x = 0:h:1
U = η(x)
t = 0.0
plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
plot!([0,1],[U[1],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = TRBDF2_heat(U)
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot(x, U, xaxis = [0,1], yaxis = [-1,2],lw=3,label = @sprintf("u(x,t), t = %1.2f",t))
        plot!([0,1],[U[1],U[end]], label = "BCs", seriestype = :scatter) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"heat_TRBDF2_noflux.gif")

η = (x,y) -> exp.(-20*(x .-1/2).^2 -20*(y .-1/2).^2)

x = 0:h:1 |> Array
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x));

function plot_heat(U,x,y,t,cl,width=800)
    l = @layout [a b]
    p1 = surface(x, y, U[end:-1:1,:], zaxis = [cl[1],cl[2]], clims= cl, aspectratio = .6, xlabel = L"x", ylabel = L"y", zlabel = L"u(x,y,t)")
    p2 = contour(x, y, U[end:-1:1,:], clims=cl, xaxis = [0,1], yaxis = [0,1], fill = true, aspectratio = 1, xlabel = L"x", ylabel = L"y")
    plot(p1, p2, layout = 2, size = (width, 7*width/10), title = @sprintf("t = %1.4f",t))
end

plot()
anim = Animation()
n = convert(Int64,ceil(T/k))

U = η(X,Y)
t = 0.0
plot_heat(U,x,y,t,(-0.5,1.5)) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = TRBDF2_heat(U) # diffuse in y
    U = TRBDF2_heat(U')' |> Array # diffuse in x
    # Since the matrix $B$ is the same if we reverse the order
    # of the grid points, we don't have to worry about what the
    # transpose does to the order.
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(U,x,y,t,(-0.5,1.5)) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"heat_TRBDF2_noflux_2D.gif")

η = (x,y) -> exp.(-20*(x .-1/2).^2 -20*(y .-1/2).^2)
F = (x,y) -> 16*sin.(4*pi*x).*sin.(4*pi*y)

function TRBDF2_heat(U,F1,F2)
    Ustar = B1\(B2*U + F1)
    return B3\((4/3)*Ustar-(1/3)*U + F2)
end

plot()
anim = Animation()
T = 2
n = convert(Int64,ceil(T/k))

U = η(X,Y)
Fh = F(X,Y)
F1 = (k/4)*Fh # γ = 1/2
F2 = (k/6)*Fh
t = 0.0
plot_heat(U,x,y,t,(-0.5,1.5)) |> IJulia.display
frame(anim)

fr = 100 #frames/unit time
tb = convert(Int64,ceil(n/(fr*T)))
for i = 2:n+1
    t += k
    U = TRBDF2_heat(U,F1,F2)
    U = TRBDF2_heat(U',F1',F2')' |> Array # Could use Strang splitting for second order in time
    if mod(i-1,tb) ≈ 0.0
        IJulia.clear_output(true)
        plot_heat(U,x,y,t,(-0.5,1.5)) |> IJulia.display
        frame(anim)
    end
end
gif(anim,"heat_TRBDF2_noflux_forcing_2D.gif")
