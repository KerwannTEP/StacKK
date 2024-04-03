# include("Args.jl")
# include("Constants.jl")
# include("ReducedVariables.jl")

using PolynomialRoots

function KroneckerDelta(i::Int64, j::Int64)
    if (i==j)
        return 1
    else
        return 0
    end
end

function _ri(i::Int64)

    return -KroneckerDelta(1,i)-KroneckerDelta(2,i)+6*KroneckerDelta(3,i)-12*KroneckerDelta(4,i)+12*KroneckerDelta(5,i)
end

function _sigmai(i::Int64)

    return 1+KroneckerDelta(1,i)+2*KroneckerDelta(2,i)+3*KroneckerDelta(3,i)+4*(KroneckerDelta(4,i)+KroneckerDelta(5,i))
end

function w_eps(v::Float64, eps::Float64, z::Float64)

    return v^2+2*eps*sqrt(z)*v+1.0
end 

function q_eps(v::Float64, eps::Float64, z::Float64, tE::Float64)

    return v^4 + 2.0*(eps*sqrt(z)+2*ta*tE)*v^3+2.0*v^2+2.0*(eps*sqrt(z)-2.0*ta*tE)*v + 1.0
end

function alpha_w(eps::Float64, z::Float64)

    solutions = roots([1.0, 2.0*eps*sqrt(z), 1.0])

    return solutions
end

function beta_q(eps::Float64, z::Float64, tE::Float64)

    solutions = roots([1.0, 2.0*(eps*sqrt(z)-2.0*ta*tE), 2.0, 2.0*(eps*sqrt(z)+2.0*ta*tE), 1.0])

    return solutions 
end

function x_eps(tE::Float64, eps::Float64, z::Float64, t::Float64)

    return 2.0*ta*tE*t*sqrt(1.0-t^2)/(1.0+eps*sqrt(z)*t)
end


# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function F_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=100)

    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    z = 2*h*tE*tLz^2

    sum = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum += num/den 
    end 

    sum *= f0*tE^(5/2)*1.0/nbK

    return sum 
end 

function F(E::Float64, Lz::Float64, nbK::Int64=100)

    return F_eps(-1.0,E,Lz,nbK) + F_eps(1.0,E,Lz,nbK)

end



function F_actions(Ju::Float64, Jv::Float64, Lz::Float64, nbK::Int64=100)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    return F(E,Lz,nbK)

end

# # atan Ju, atan Jv, atan Lz
# function get_DF_actions(nbJ::Int64=100, nbK::Int64=100)

#     tab_F_JuJvLz = zeros(Float64, nbJ, nbJ, nbJ)

#     for iu=1:nbJ 
#         atanJu = pi/2*(iu-1)/(nbJ -1)
#         Ju = tan(atanJu)
#         for iv=1:nbJ 
#             atanJv = pi/2*(iv-1)/(nbJ -1)
#             Jv = tan(atanJv)
#             for iz=1:nbJ 
#                 atanLz = -pi/2+pi*(iz-1)/(nbJ -1)
#                 Lz = tan(atanLz)
#                 println((Ju,Jv,Lz))
#                 # E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)
#                 if (iu != nbJ && iv != nbJ && iz != 1 && iz != nbJ)
#                     tab_F_JuJvLz[iu,iv,iz] = F_actions(Ju,Jv,Lz,nbK)
#                 else
#                     tab_F_JuJvLz[iu,iv,iz] = 0.0
#                 end
#             end
#         end

#     end

#     return tab_F_JuJvLz
# end

# const nbJsampling = 100

# const tab_DF = get_DF_actions(nbJsampling)

# function itp_DF_actions(Ju::Float64, Jv::Float64, Lz::Float64, nbJ::Int64=nbJsampling)

#     atanJu = atan(Ju)
#     atanJv = atan(Jv)
#     atanLz = atan(Lz)

#     iu = 1 + floor(Int64,atanJu/(pi/2.0/(nbJ -1)))
#     iv = 1 + floor(Int64,atanJv/(pi/2.0/(nbJ -1)))
#     iz = 1 + floor(Int64,(atanLz+pi/2.0)/(pi/(nbJ -1)))

#     # Bilinear/Trilinear

#     return tab_D[iu,iv,iz]

# end



# Derivatives

# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdE_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=100)


    # dF/dtE
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    z = 2*h*tE*tLz^2

    # Derivative on tE^2.5
    sum1 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum1 += num/den 
    end 

    sum1 *= 2.5*f0*tE^(3/2)*1.0/nbK

    dzdtE = 2*h*tLz^2

    # Derivative within the integral 
    sum2 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        dxedtE = 2*ta*sqrt(1-t^2)/(1+eps*sqrt(z)*t) + 2*ta*tE*t*sqrt(1-t^2)*(-eps*t*0.5*dzdtE/sqrt(z))/(1+eps*sqrt(z)*t)^2



        # num = (1 - t^2) * [ (3 + 4 xe - xe^2)(1 - xe)(1-t^2) + 12 t^2 ]
        # den = (1 - 2 ta tE t sqrt(1-t^2) + eps sqrt(z) t)^5

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        dnum = (1-t^2)*((4*dxedtE-2*xe*dxedtE)*(1-xe)*(1-t^2) + (3+4*xe-xe^2)*(-dxedtE)*(1-t^2))
        dden = 5 * (-2*ta*t*sqrt(1-t^2) + eps*t*(0.5*dzdtE/sqrt(z)) ) * (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^4


        sum2 += (dnum*den-num*dden)/den^2 
    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK


    # dF/dE = dF/dtE dtE/dE
    # dtE/dE = -(a+c)/(G*M)

    sum = sum1 + sum2
    sum *= -(a+c)/(G*M)

    return sum 
end 


function _dFdE(E::Float64, Lz::Float64, nbK::Int64=100)

    return _dFdE_eps(-1.0,E,Lz,nbK) + _dFdE_eps(1.0,E,Lz,nbK)

end



# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdLz_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=100)


    # dF/dtLz
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    z = 2*h*tE*tLz^2

   
    dzdtLz = 2*h*tE*2*tLz

    # Derivative within the integral 
    sum2 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        dxedtLz =  2*ta*tE*t*sqrt(1-t^2)*(-eps*t*0.5*dzdtLz/sqrt(z))/(1+eps*sqrt(z)*t)^2



        # num = (1 - t^2) * [ (3 + 4 xe - xe^2)(1 - xe)(1-t^2) + 12 t^2 ]
        # den = (1 - 2 ta tE t sqrt(1-t^2) + eps sqrt(z) t)^5

        num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        dnum = (1-t^2)*((4*dxedtLz-2*xe*dxedtLz)*(1-xe)*(1-t^2) + (3+4*xe-xe^2)*(-dxedtLz)*(1-t^2))
        dden = 5 * ( eps*t*(0.5*dzdtLz/sqrt(z)) ) * (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^4


        sum2 += (dnum*den-num*dden)/den^2 
    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK


    # dF/dLz = dF/dtLz dtLz/dLz
    # dtLz/dELz = 1/sqrt((a+c)*G*M)


    sum2 *= 1/sqrt((a+c)*G*M)

    return sum2
end 

function _dFdLz(E::Float64, Lz::Float64, nbK::Int64=100)

    return _dFdLz_eps(-1.0,E,Lz,nbK) + _dFdLz_eps(1.0,E,Lz,nbK)

end


















function test_dFdE(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    E_p = E + eps 
    E_m = E - eps 

    dFdE_num = (F(E_p,Lz,nbK) - F(E_m,Lz,nbK))/(2*eps)
    dFdE_th = _dFdE(E,Lz,nbK)

    println("Num = ",dFdE_num)
    println("Th  = ",dFdE_th)


end


function test_dFdLz(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    Lz_p = Lz + eps 
    Lz_m = Lz - eps 

    dFdLz_num = (F(E,Lz_p,nbK) - F(E,Lz_m,nbK))/(2*eps)
    dFdLz_th = _dFdLz(E,Lz,nbK)

    println("Num = ",dFdLz_num)
    println("Th  = ",dFdLz_th)


end







function integrate_F_actions(nbJu::Int64, nbJv::Int64, nbLz::Int64)

    # use Ji = tan(J)
    # dJi = (1+tan^2(J)) dJ 

    nbgrid = nbJu*nbJv*nbLz 

    index = 1
    tabJuJvJLz = zeros(Float64, 3, nbgrid)

    for i=1:nbJu 

        Ju = pi/2 * (i-0.5)/nbJu

        for j=1:nbJv 

            Jv = pi/2 * (j-0.5)/nbJv

            for k=1:nbLz

                Lz = pi/2 * (k-0.5)/nbLz
                tabJuJvJLz[1,index], tabJuJvJLz[2,index], tabJuJvJLz[3,index] = Ju, Jv, Lz
            
                index += 1
            end

        end

    end



    sum = Threads.Atomic{Float64}(0.0)

    Threads.@threads for igrid=1:nbgrid 
        Ju, Jv, Lz = tabJuJvJLz[1,igrid], tabJuJvJLz[2,igrid], tabJuJvJLz[3,igrid]
        jacJu = 1 + tan(Ju)^2
        jacJv = 1 + tan(Jv)^2
        jacLz = 1 + tan(Lz)^2
        FJ = F_actions(tan(Ju),tan(Jv),tan(Lz))

        Threads.atomic_add!(sum, jacLz*pi/2/nbLz*FJ*jacJv*pi/2/nbJv *jacJu*pi/2/nbJu  )
    end

    sum[] *= 2.0*(2*pi)^3 

    return sum[]


end




