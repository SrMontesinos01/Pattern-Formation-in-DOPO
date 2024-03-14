module stability

using FFTW
import Plots

σ(k::Float64, E::Float64) = -1 +sqrt(E^2 - (Δ1 + 2*k^2)^2)
E_f(k::Float64) = sqrt(1 + (Δ1 + 2k^2)^2)
Δ1 = -0.18
kc = sqrt(-Δ1/2.0)
λc = 2*π/ kc
L = 4*λc

E = 0.999
maxval = sqrt((E -Δ1)/2.0)
# -------------------------- Dispersion Plot -------------------------- #
Δx1 = 1.32
Δx2 = 3.0
if true
    xArr1 = collect(0:Δx1:L)
    kArr1 = 2π * fftshift(fftfreq(length(xArr1), length(xArr1)/L))
    kArr1 = filter(x -> abs(x) < maxval, kArr1)
    xArr2 = collect(0:Δx2:L)
    kArr2 = 2π *fftshift(fftfreq(length(xArr2), length(xArr2)/L))
    kArr2 = filter(x -> abs(x) < maxval, kArr2)

    E = 0.999
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E")
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, label = "",markercolor=:red)
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, label = "",markercolor=:green)

    E = 1.00
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E", color =:orange)
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, markercolor=:green, label = "")

    E = 1.001
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E")
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4,markercolor=:green, label = "")

    E = 1.01
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E", color =:purple)
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, markercolor=:green, label = "")
    # Plots.scatter!(title = "\$Δ_1 = -0.18\$, \$E = $E\$")
    Plots.scatter!(title = "Dispersion Relation \$λ_1(k)\$")
    Plots.scatter!(ylabel = "\$λ_1(k)\$")
    # root = joinpath(@__DIR__, "Images", "dispersion.png")

    E = 1.1
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E", color=:yellow)
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label="\$L = 4λ_c\$")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, markercolor=:green, label="\$L =2.5λ_c\$")
    Plots.hline!([0], linestyle=:dash, color=:red, label="\$λ_1(k) = 0\$")
end

Plots.xlims!(0.0, 1)  
Plots.ylims!(-0.005, 0.011)
Plots.savefig("dispersion.png")
# -------------------------- Dispersion Plot v2 -------------------------- #
L1 = 4*λc
L2 = 2.5*λc
n = 64

if true
    Δx1 = L1 / (n - 1)
    Δx2 = L2 / (n - 1)

    xArr1 = collect(0:Δx1:L)
    kArr1 = 2π * fftshift(fftfreq(length(xArr1), length(xArr1)/L1))
    kArr1 = filter(x -> abs(x) < maxval, kArr1)
    xArr2 = collect(0:Δx2:L)
    kArr2 = 2π *fftshift(fftfreq(length(xArr2), length(xArr2)/L2))
    kArr2 = filter(x -> abs(x) < maxval, kArr2)

    E = 0.999
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E")
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, label = "",markercolor=:red)
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, label = "",markercolor=:green)

    E = 1.00
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E", color =:orange)
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, markercolor=:green, label = "")

    E = 1.001
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E")
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4,markercolor=:green, label = "")

    E = 1.005
    k_range = 0:0.01:0.65
    y_range = σ.(k_range, E)
    puntos_x1 = kArr1
    puntos_y1 = σ.(puntos_x1, E)
    puntos_x2 = kArr2
    puntos_y2 = σ.(puntos_x2, E)

    Plots.plot!(k_range, y_range, linewidth=2, xlabel = "k", label = "\$E = \$ $E", color =:purple)
    Plots.hline!([0], linestyle=:dash, color=:red, label="\$λ_+(k) = 0\$")
    Plots.scatter!(puntos_x1, puntos_y1, markersize=4, markercolor=:red, label = "\$L =4λ_c\$")
    Plots.scatter!(puntos_x2, puntos_y2, markersize=4, markercolor=:green, label = "\$L =2.5λ_c\$")
    # Plots.scatter!(title = "\$Δ_1 = -0.18\$, \$E = $E\$")
    Plots.scatter!(title = "Dispersion Relation \$λ_+(k)\$")
    Plots.scatter!(ylabel = "\$λ_+(k)\$")
    
    
end

Plots.xlims!(0.1, 0.5)  
Plots.ylims!(-0.005, 0.006)
root = joinpath(@__DIR__, "Images", "dispersion.png")
Plots.savefig(root)
# -------------------------- Marginal stability Plot -------------------------- #
k_range = 0:0.001:0.9
y_range = E_f.(k_range)
Plots.plot(k_range , y_range, label="f(x)", linewidth=2)
Plots.hline!([1], linestyle=:dash, color=:green, label="E = 1.0")
Plots.hline!([1.0001], linestyle=:dash, color=:orange, label="E = 1.0001")
Plots.hline!([1.0005], linestyle=:dash, color=:red, label="E = 1.005")
Plots.hline!([1.05], linestyle=:dash, color=:red, label="E = 1.05")

Δx = 1.0
xArr = collect(0:Δx:L)
kArr = fftshift(fftfreq(length(xArr)))*2π
puntos_x = kArr # Points that are taken in the discretization
puntos_y = E_f.(puntos_x)
Plots.scatter!(puntos_x, puntos_y, markersize=2, markercolor=:purple, label="Puntos específicos")

Plots.xlims!(0.0, 0.9)  
Plots.ylims!(0.98, 1.3) 
end