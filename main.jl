module mainV2

include("auxFunc.jl")

using DelimitedFiles, FFTW, .auxFunc
using GLMakie # Use the activate DynSyst_Project

# ----------------------------------------------------------- #
# ------------------   System Parameters   ------------------ #
# ----------------------------------------------------------- #
Δ0, Δ1, c = 0.0, -0.18, 10^(-4) # g/(√aγ)

kc = sqrt(-Δ1/2.0)
λc = 2*π/ kc
L = 4*λc

z0, z1, = 1 +im*Δ0, 1 +im*Δ1
a0, a1 = 1.0, 2.0
b = sqrt(2)*c
b = 0
# ----------------------------------------------------------- #
# ------------------ Simulation Parameters ------------------ #
# ----------------------------------------------------------- #
Δx = 1.32
Δt = 0.01
nwrt = 20 # usual value = 200
nstp = 10^5
totalTime = nwrt*nstp*Δt
n = 64 # Constant number of points of the simulation. 

parameters = collect([Δx, Δt, z0, z1, a0, a1, b, nwrt, nstp])
xArr = collect(0:Δx:L)
tArr = collect(0:Δt:totalTime)

Δx = 0.3
2*a0*Δt/ (Δx)^2 # Courant criterions 
2*a1*Δt/ (Δx)^2 # Courant criterions 
# ----------------------------------------------------------- #
# ------------------       Simulation      ------------------ #
# ----------------------------------------------------------- #
listE = [0.999+ im*0.0, 1.0 + im*0.0, 1.001 + im*0.0, 1.01 + im*0.0, 1.1 + im*0.0]
listL = [3*λc,  4*λc, 6*λc]

listE = [0.999 + im*0.0, 1.0 + im*0.0, 1.001 + im*0.0, 1.005 + im*0.0,]
listL = [2.5*λc, 4*λc]
listL = [4.0*λc]

# listΔx = listL / (n -1)
# Δx = round(4.0*λc / (n-1); digits = 2)
# xArr = collect(0:listΔx[1]:listL[1])

# It takes ≈3h30mins for 4×5 values of E×Δx with nwrt = 200, nstp = 10^5
@time begin
    for L in listL
        Δx = round(L / (n - 1); digits = 2)
        print("----------------- L = $L ----------------- \n")
        parameters = collect([Δx, Δt, z0, z1, a0, a1, b, nwrt, nstp])
        for E in listE
            print("E = $E,")

            
            xArr = collect(0:Δx:L)
            A0 = (E / z0)*complex.(ones(length(xArr))) # Initial Condition
            A1 = complex.((10^(-5))*(10*cos.(kc*xArr) + randn(length(xArr)))) # Initial Condition
            A1 = complex.(zeros(length(xArr)))
            t = 0.0 # Initial Condition

            integrationSPDE(A0, A1, t, E, parameters)
        end 
        print("\n") 
    end
end

# ----------------------------------------------------------- #
# ------------------        Ploting        ------------------ #
# ----------------------------------------------------------- #
function plotHeatMap(E, Δx)
    root = joinpath(@__DIR__, "Data", "A1 $Δx $E.txt")
    Data = readdlm(root, '\t', ComplexF64)
    
    matrix = real.(Data[:,2:end])
    matrixFFT = abs.(mapslices(fft, Data[:,2:end], dims=2)).^2
    matrixFFT = log.(mapslices(fftshift, matrixFFT, dims=2))
    
    tArr = abs.(Data[:,1])
    # print(argmax(matrixFFT, dims=2))
    Δx = real(Δx)
    xArr = collect(0:Δx:L)
    # kArr = fftshift(fftfreq(length(xArr), Δx))
    kArr = 2π * fftshift(fftfreq(length(xArr), length(xArr)/L))

    fig = Figure(figure_padding = 15)
    ax1 = Axis(fig[1,1]; xlabel = L"x", title = L"α_{1R}(x,t)", ylabel = L"\text{Time} (t)", yticks =[false])
    ax2 = Axis(fig[1,2]; xlabel = L"k", title = L"|α_1(k,t)|^2", yticks =[false])
    ax1.xlabelsize = 20
    ax1.ylabelsize = 20
    ax2.xlabelsize = 20
    ax2.titlesize = 30
    ax1.titlesize = 30
    hm1 = heatmap!(ax1, xArr, tArr, matrix', colormap =:oslo)
    hm2 = heatmap!(ax2, kArr, tArr, matrixFFT', colormap =:oslo)

    # hidexdecorations!(ax2, ticks=false)
    Colorbar(fig[:, end+1], hm1, label = "NF")
    # Colorbar(fig[:, end+1], hm2, label = "FF")

    Etitle = real(E)
    tit = L"\text{Intensity of the pump field } E = %$Etitle \text{ }(L = 4λ_c)"
    fig[0, :] = Label(fig, tit , fontsize = 30)

    root = joinpath(@__DIR__, "Images", "plot E_$Etitle x_$Δx v2.png")
    save(root, fig, px_per_unit = 2)
end

listE = [0.999 + im*0.0, 1.0 + im*0.0, 1.001 + im*0.0, 1.01 + im*0.0, 1.1 + im*0.0]
listL = [4*λc]

listE = [0.999 + im*0.0, 1.0 + im*0.0, 1.001 + im*0.0, 1.005 + im*0.0,]
listL = [2.5*λc, 4*λc]
listL = [4*λc]

# Loop for creating and saving the figures
for L in listL
    for E in listE
        Δx = round(L / (n - 1); digits = 2)
        plotHeatMap(E, complex(Δx))
    end
end

Δx = 1.32
xArr = collect(0:Δx:L)
# kArr = fftshift(fftfreq(length(xArr), Δx))
kArr = 2π *fftfreq(length(xArr)) # [0:n÷2-1; -n÷2:-1]  * (fs/n) *(2π), fs = 1 by default
kArr = fftshift(kArr)
kArr[47]
kArr[39]

"""
E = listE[1]
Δx= complex(listΔx[1])
root = joinpath(@__DIR__, "Data", "A1 $Δx $E.txt")
Data = readdlm(root, '\t', ComplexF64)

matrix = real.(Data[:,2:end])
matrixFFT = abs.(mapslices(fft, Data[:,2:end], dims=2)).^2
matrixFFT = real.(abs.(mapslices(fft, Data[:,2:end], dims=2)))
matrixFFT = log.(mapslices(fftshift, matrixFFT, dims=2))

Δx= complex(listΔx[1])
Δx = real(Δx)
xArr = collect(0:Δx:L)
kArr = 2π * fftshift(fftfreq(length(xArr)))
kArr[47]
kArr[39]
"""
# ------------------------------------------------------------ #
function fftfreqHomeMade(L::Float64, Δx::Float64)
    xArr = 0:Δx:L
    N = length(xArr)
    kArr = (2π/L) * [0:N÷2-1; -N÷2:-1]
    return kArr
end

function fftHomeMade(v::Vector{ComplexF64})
    xArr = 0:Δx:L
    N = length(xArr)
    kArr = fftfreqHomeMade(L, Δx)

    vk = complex(zeros(N))
    for l in 1:N
        k = kArr[l]
        for j in 1:N
            x = xArr[j]
            vk[l] += v[j] * exp(im*k*x)
        end
    end

    return vk
end

# This one uses the "HomeMade" functions
function plotHeatMap_v2(E, Δx)
    root = joinpath(@__DIR__, "Data", "A1 $Δx $E.txt")
    Data = readdlm(root, '\t', ComplexF64)
    
    matrix = real.(Data[:,2:end])
    # matrixFFT = abs.(mapslices(fft, Data[:,2:end], dims=2)).^2
    # matrixFFT = log.(mapslices(fftshift, matrixFFT, dims=2))

    matrixFFT = abs.(mapslices(fftHomeMade, Data[:,2:end], dims=2)).^2
    matrixFFT = log.(matrixFFT)
    tArr = abs.(Data[:,1])
    
    Δx = real(Δx)
    xArr = collect(0:Δx:L)
    # kArr = fftshift(fftfreq(length(xArr), 1 /Δx))
    kArr = fftfreqHomeMade(L, Δx)

    fig = Figure(figure_padding = 15)
    ax1 = Axis(fig[1,1]; xlabel = L"x", title = L"α_{1R}(x,t)", ylabel = L"\text{Time} (t)", yticks =[false])
    ax2 = Axis(fig[1,2]; xlabel = L"k", title = L"|α_1(k,t)|^2", yticks =[false])
    ax1.xlabelsize = 20
    ax1.ylabelsize = 20
    ax2.xlabelsize = 20
    ax2.titlesize = 30
    ax1.titlesize = 30
    hm1 = heatmap!(ax1, xArr, tArr, matrix', colormap =:oslo)
    hm2 = heatmap!(ax2, kArr, tArr, matrixFFT', colormap =:oslo)

    # hidexdecorations!(ax2, ticks=false)
    # Colorbar(fig[:, end+1], hm1, label = "NF")
    # Colorbar(fig[:, end+1], hm2, label = "FF")

    Etitle = real(E)
    tit = L"\text{Intensity of the pump field } E = %$Etitle \text{ }(Δx = %$Δx)"
    fig[0, :] = Label(fig, tit , fontsize = 30)

    root = joinpath(@__DIR__, "Images", "plot E_$Etitle x_$Δx V2.png")
    save(root, fig, px_per_unit = 2)
end

# Loop for creating and saving the figures using HOMEMADE
for Δx in complex.(listΔx)
    for E in listE
        plotHeatMap_v2(E, Δx)
    end
end

k = (2π/L) * [0:N÷2-1; -N÷2:-1]

N = 100  # ajusta según sea necesario
Δx = 0.1
L = N * Δx
x = 0:Δx:L
x = 2π*(0:N-1)/N


Δx = 2π / N
X = fft(x)
k = fftfreq(N, abs(x[2] - x[1]))
Δk = 2π / L

Δk, k[3] - k[2]

n = N
fs = Δx
[0:n÷2-1; -n÷2:-1]  * fs/n 


σ(k::Float64, E::Float64) = -1 + sqrt(E^2 -(Δ1 + 2*k^2))
Δx = 1.32
xArr = collect(0:Δx:L)
kArr = fftshift(fftfreq(length(xArr), Δx))
kArr[33:end]
x_range = 0:0.01:0.65
y_range = σ.(x_range, 0.999)

fig = lines(x_range, y_range, color=:blue, linewidth=2, label="f(x)")

x_range = kArr[33:end]
y_range = σ.(kArr[33:end], 0.999)
scatter!(fig, x_range, y_range, markersize=8, markercolor=:red, label="Puntos específicos")



"""
x = tArr
y1 = abs.(Data[:,2]).^2
y2 = abs.(Data[:,60]).^2
y3 = abs.(Data[:,20]).^2
# Crear el gráfico de dispersión
fig = Figure()
ax = Axis(fig[1,1]; xlabel = "time", ylabel = "variable")
lines!(ax, x, y1)
lines!(ax, x, y2)
lines!(ax, x, y3)
"""

# ------------------- testing stuff ------------------- #
using GLMakie

# Define tu función f(x)
f(x) = sin(x)

# Genera un rango de valores x
x_range = 0:0.1:2π

# Calcula los valores de la función f(x)
y_values = f.(x_range)

# Crea el gráfico con la línea de la función
fig, ax = lines(x_range, y_values, color=:blue, linewidth=2, label="f(x)")

# Define los puntos específicos
puntos_x = [π/2, 3π/2]
puntos_y = f.(puntos_x)

# Añade marcadores para los puntos específicos
scatter!(ax, puntos_x, puntos_y, markersize=8, markercolor=:red, label="Puntos específicos")

# Muestra el gráfico
display(fig)

Δq = 2π / (N*Δx)
Δq = 2π / L
qmax = 2π/ Δx
qArr = collect(0:Δq:qmax)
fftshift(qArr)

L 
Δx
kArr, vk = fftHomeMade(L, Δx, complex(xArr))
kArr
kArr[2] - kArr[1]
vk


a = 0.5 - 1im
b = 5 + 5im
c = 3im
d = 4
typeof(a)
abs(a)
real(a)
Int64(real(5.0 + im*0))
convert(Int, real(5.0 + im*0))

randn(2)

v = collect([0+4im,1-2im,2,1])
w = fft(v)

a = conj(a)
v1 = collect([a,b,c,d])
v2 =  collect([a,b,c,d])
v1.*a
length(v1)
# f0.(v1, a)
# f1.(v1, v2)
# wNoise.(1, v1)

v1 + v2
shift(v1, 1)
v1
# sp = ShiftedArrays.circshift(v1, 1)
#sm = ShiftedArrays.circshift(v1, -1)

v = complex(collect([1,2,3,4,5]))
w = complex(collect([1,2,3,4,5]))
(1 + im)*v + 0.0*w
z = 50.0 + 6*im
Int64(real(z))
end