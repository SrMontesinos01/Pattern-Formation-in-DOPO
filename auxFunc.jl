module auxFunc
using  ShiftedArrays, DelimitedFiles

export integrationSPDE

function g1(A_0::ComplexF64)
    #print(A_0)
    num1 = -imag(A_0) / (2*sqrt(2 + real(A_0)))
    num2 = im*0.5*sqrt(2 + real(A_0))

    return num1 + num2
end

function g2(A_0::ComplexF64)
    if (abs(A_0) < 2)
        num1 = 1 - (abs(A_0)^2)/4.0
        num2 = 2 + real(A_0)

        return sqrt(num1/num2)
    else
        print("Error abs(A_0) > 2")
        return 0
    end
end

function f0(A_1::ComplexF64, E::ComplexF64)
    return E - 0.5*A_1*A_1
end

function f1(A_0::ComplexF64, A_1::ComplexF64)
    return A_0*conj(A_1)
end

function wNoise(index::Int64, A_0::ComplexF64)
    if index == 0
        return randn()
    elseif index == 1
        return g1(A_0)*randn() + g2(A_0)*randn()
    else
        print("Error: Index must be 0 or 1")
        return 0
    end
end

function heunMethodSPDE(A0::Vector{ComplexF64}, A1::Vector{ComplexF64}, t::Float64, E::ComplexF64, params::Vector{})
    diff, D0, D1, Δt, Δt2, z0, z1, b, N = params

    # Shifting Arrays for periodic boundary conditions
    A0_p = ShiftedArrays.circshift(A0, 1)
    A0_m = ShiftedArrays.circshift(A0, -1)
    A1_p = ShiftedArrays.circshift(A1, 1)
    A1_m = ShiftedArrays.circshift(A1, -1)
    
    N = Int64(real(N)) # This should be the number of points we take in the discretization
    u1, u2, u3 = randn(N), randn(N), randn(N) # Random numbers for the corresponding step

    # Noise term for the "predictor" (Implicit euler)
    w0 = b*diff*u1
    w1 = b*diff*(g1.(A0).*u2+ g2.(A0).*u3)

    # Computing the "predictor" (Implicit euler)
    λ0 = Δt*(-z0*A0 + D0*(A0_p - 2*A0 + A0_m) + f0.(A1, E)) + w0
    λ1 = Δt*(-z1*A1 + D1*(A1_p - 2*A1 + A1_m) + f1.(A0, A1)) + w1
    A0_1 = A0 + λ0
    A1_1 = A1 + λ1

    # Shifting Arrays for periodic boundary conditions
    A0_1p = ShiftedArrays.circshift(A0_1, 1)
    A0_1m = ShiftedArrays.circshift(A0_1, -1)
    A1_1p = ShiftedArrays.circshift(A1_1, 1)
    A1_1m = ShiftedArrays.circshift(A1_1, -1)

    # Noise term for the last step
    w1 = b*diff*(g1.(A0_1).*u2+ g2.(A0_1).*u3)

    # Compute the final new values for the fields
    A0_new = A0 + 0.5*λ0 + Δt2*(-z0*A0_1 .+ D0*(A0_1p - 2*A0_1 + A0_1m) + f0.(A1_1, E)) .+ 0.5*w0
    A1_new = A1 + 0.5*λ1 + Δt2*(-z1*A1_1 .+ D1*(A1_1p - 2*A1_1 + A1_1m) + f1.(A0_1, A1_1)) .+ 0.5*w1
    t = t + Δt

    return A0_new, A1_new, real(t)
end

function writeFile(A::Vector{ComplexF64}, t::Float64, f::IOStream)
    T = t + im*0.0
    write(f, "$T")
    for num in A
        write(f, "\t$num")
    end
    write(f, "\n")
end

function integrationSPDE(A0::Vector{ComplexF64}, A1::Vector{ComplexF64}, t::Float64, E::ComplexF64, params::Vector{})
    Δx, Δt, z0, z1, a0, a1, b, nwrt, nstp = params
    diff = sqrt(Δt/Δx)
    D0 = (im*(a0)/((Δx)^2))
    D1 = (im*(a1)/((Δx)^2))
    Δt2 = Δt/2.0
    N = length(A0)
    p = collect([diff, D0, D1, Δt, Δt2, z0, z1, b, N])
    # print(p)

    root0 = joinpath(@__DIR__, "Data" ,"A0 $Δx $E.txt")
    root1 = joinpath(@__DIR__, "Data" ,"A1 $Δx $E.txt")
    file0 = open(root0, "w")
    file1 = open(root1, "w")

    nwrt = convert(Int, real(nwrt))
    nstp = convert(Int, real(nstp))
    for i in 1:nwrt
        # @show i
        if i % 10 == 0
            print("$i, ")
        end
        writeFile(A0, t, file0)
        writeFile(A1, t, file1) 
        for j in 1:nstp
            A0, A1, t = heunMethodSPDE(A0, A1, t, E, p)
        end
    end
    close(file0)
    close(file1)
end
    
end