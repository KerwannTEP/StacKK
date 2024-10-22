include("Main.jl")

function lambda0_lambda_1(Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=nbu_default)

    E, _, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv,nbu)
    u0, u1 = find_bounds_u(E,Lz,I3)

    lambda0 = lambda_from_u(u0)
    lambda1 = lambda_from_u(u1)

    return lambda0, lambda1
end


function diff_u(Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=nbu_default)

    E, _, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv,nbu)
    u0, u1 = find_bounds_u(E,Lz,I3)

    return u1-u0
end

function diff_lambda(Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=nbu_default)

    lambda0, lambda1 = lambda0_lambda_1(Ju,Jv,Lz,nbu)

    return lambda1 - lambda0
end

using Plots

function plot_E(Jv::Float64, Lz::Float64, nbu::Int64=nbu_default)

    tabJu = [i/100 for i=1:100]
    tabE = [E_Lz_I3_from_Ju_Lz_Jv(tabJu[i],Lz,Jv,nbu)[1] for i=1:100]
   
    plot(tabJu,tabE)
end

function plot_diff_u(Jv::Float64, Lz::Float64, nbu::Int64=nbu_default)

    tabJu = [i/1000 for i=1:1000]
    tabdiff = [diff_u(tabJu[i],Lz,Jv,nbu) for i=1:1000]
   
    plot(tabJu,tabdiff)
end