# include("Args.jl")
# include("Constants.jl")
# include("ReducedVariables.jl")

using PolynomialRoots

# function KroneckerDelta(i::Int64, j::Int64)
#     if (i==j)
#         return 1
#     else
#         return 0
#     end
# end

# function _ri(i::Int64)

#     return -KroneckerDelta(1,i)-KroneckerDelta(2,i)+6*KroneckerDelta(3,i)-12*KroneckerDelta(4,i)+12*KroneckerDelta(5,i)
# end

# function _sigmai(i::Int64)

#     return 1+KroneckerDelta(1,i)+2*KroneckerDelta(2,i)+3*KroneckerDelta(3,i)+4*(KroneckerDelta(4,i)+KroneckerDelta(5,i))
# end

# function w_eps(v::Float64, eps::Float64, z::Float64)

#     return v^2+2*eps*sqrt(z)*v+1.0
# end 

# function q_eps(v::Float64, eps::Float64, z::Float64, tE::Float64)

#     return v^4 + 2.0*(eps*sqrt(z)+2*ta*tE)*v^3+2.0*v^2+2.0*(eps*sqrt(z)-2.0*ta*tE)*v + 1.0
# end

# function alpha_w(eps::Float64, z::Float64)

#     solutions = roots([1.0, 2.0*eps*sqrt(z), 1.0])

#     return solutions
# end

# function beta_q(eps::Float64, z::Float64, tE::Float64)

#     solutions = roots([1.0, 2.0*(eps*sqrt(z)-2.0*ta*tE), 2.0, 2.0*(eps*sqrt(z)+2.0*ta*tE), 1.0])

#     return solutions 
# end

function x_eps(tE::Float64, eps::Float64, z::Float64, t::Float64)

    return 2.0*ta*tE*t*sqrt(1.0-t^2)/(1.0+eps*sqrt(z)*t)
end


# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function F_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)

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



function F_actions(Ju::Float64, Jv::Float64, Lz::Float64, nbK::Int64=nbK_default)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    return F(E,Lz,nbK)

end

using Plots

# Derivatives


# Mistake somewhere here for E-gradient?
# Plot derivative 
# Compare each components


# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdE_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)


    # dF/dtE
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2.0^(3/2)*pi^3*ta)
    z = 2.0*h*tE*tLz^2

    # Derivative on tE^2.5
    sum1 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        num = (1.0-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1.0-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        sum1 += num/den 
    end 

    sum1 *= 2.5*f0*tE^(3/2)*1.0/nbK

    dzdE = 2.0*h*tLz^2

    # Derivative within the integral 
    # ISSUE HERE ?

    # data = zeros(Float64,nbK)

    sum2 = 0.0
    for k=1:nbK
        
        
        # t = 1.0/nbK*(k-0.5)
        # xe = x_eps(tE,eps,z,t)

        # dxedE = 2*ta*sqrt(1.0-t^2)/(1.0+eps*sqrt(z)*t) + 2.0*ta*tE*t*sqrt(1-t^2)*(-eps*t*0.5*dzdE/sqrt(z))/(1+eps*sqrt(z)*t)^2



        
    #     ((1 - t^2)  (-2  a  t  (-1 + t^2)  (1 - 2  a  t  Sqrt[1 - t^2]  tE + 
    #      Sqrt[2]  eps  t  Sqrt[
    #       h Lz^2 tE])  (2  eps^3  h^2  Lz^4  t^3  Sqrt[2 - 2 t^2]
    #         tE^2 + 
    #      8  eps^2  t^2  (h Lz^2 tE)^(
    #       3/2)  (Sqrt[1 - t^2] + 5  a  t  (-1 + t^2)  tE) + 
    #      Sqrt[2]  eps  h  Lz^2  t  tE  (5  Sqrt[1 - t^2] - 
    #         12  a  t  (-1 + t^2)  tE  (-5 + 
    #            a  t  Sqrt[1 - t^2]  tE)) + 
    #      2  Sqrt[h Lz^2 tE]  (Sqrt[1 - t^2] - 
    #         4  a  t  (-1 + t^2)  tE  (-5 + 
    #            3  a  t  Sqrt[1 - t^2]  tE))) - 
    #   5  t  (1 + 
    #      Sqrt[2]  eps  t  Sqrt[h Lz^2 tE])  (Sqrt[2]  eps  h  Lz^2 - 
    #      4  a  Sqrt[1 - t^2]  Sqrt[
    #       h Lz^2 tE])  (12  t^2  (1 + 
    #         Sqrt[2] eps t Sqrt[h Lz^2 tE])^3 - (1 - t^2)  (1 - 
    #         2  a  t  Sqrt[1 - t^2]  tE + 
    #         Sqrt[2]  eps  t  Sqrt[
    #          h Lz^2 tE])  (-4  a^2  t^2  (-1 + t^2)  tE^2 - 
    #         8  a  t  Sqrt[1 - t^2]
    #            tE  (1 + Sqrt[2]  eps  t  Sqrt[h Lz^2 tE]) - 
    #         3  (1 + Sqrt[2] eps t Sqrt[h Lz^2 tE])^2))))/(2  Sqrt[
    # h Lz^2 tE]  (1 + Sqrt[2] eps t Sqrt[h Lz^2 tE])^4  (1 - 
    #   2 a t Sqrt[1 - t^2] tE + Sqrt[2] eps t Sqrt[h Lz^2 tE])^6);





        # num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        # den = (1-2*a*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

        # dnumdE = (4*dxedE-2*dxedE*xe)*(1-xe)*(1-t^2)+ (3+4*xe-xe^2)*(-dxedE)*(1-t^2)
        # ddendE = 5*(1-2*a*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^4 * (-2*a*t*sqrt(1-t^2) + eps*0.5*dzdE/sqrt(z)*t)



        # sum2 += (1.0-t^2)*(dnumdE*den-num*ddendE)/den^2 


        t = 1.0/nbK*(k-0.5)

        z = 2*h*tE*tLz^2
        xe = 2*ta*tE*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t)

        dzdE = 2*h*tLz^2
        dxedE = 2*ta*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t) - 2*ta*tE*t*sqrt(1-t^2)*(eps*0.5*dzdE/sqrt(z)*t)/(1+eps*sqrt(z)*t)^2


        num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        den = (1-2*ta*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

        dnumdE = (4*dxedE-2*dxedE*xe)*(1-xe)*(1-t^2)+  (3+4*xe-xe^2)*(-dxedE)*(1-t^2) 
        ddendE = 5*(1-2*ta*tE*t*sqrt(1-t^2) + eps*sqrt(z)*t)^4 * (-2*ta*t*sqrt(1-t^2) + eps*0.5*dzdE/sqrt(z)*t)

        sum2 += (1-t^2)*(dnumdE*den-num*ddendE)/den^2



        # data[k] = (1.0-t^2)*(dnumdE*den-num*ddendE)/den^2 
    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK

    # plot ?


    # dF/dE = dF/dtE dtE/dE
    # dtE/dE = -(a+c)/(G*M)

    sum = sum1 + sum2
    sum *= -(a+c)/(G*M)


    # println((-(a+c)/(G*M)*sum2))

    # pt = plot(data)

    # savefig(pt,"test.png")





    # # num 

    # dE = 0.000000001

  
    
    # tE = _tE(E+dE)
    # tLz = _tLz(Lz)
    # f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    # z = 2*h*(tE)*tLz^2

    # sump = 0.0
    # for k=1:nbK
    #     t = 1.0/nbK*(k-0.5)
    #     xe = x_eps(tE,eps,z,t)

    #     num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
    #     den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

    #     sump += num/den 
    # end 

    # tE = _tE(E-dE)
    # tLz = _tLz(Lz)
    # f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2^(3/2)*pi^3*ta)
    # z = 2*h*(tE)*tLz^2

    # summ = 0.0
    # for k=1:nbK
    #     t = 1.0/nbK*(k-0.5)
    #     xe = x_eps(tE,eps,z,t)

    #     num = (1-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1-t^2)+12.0*t^2)
    #     den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

    #     summ += num/den 
    # end 


    # sum = (sump-summ)/(2.0*dE)

    # sum *= f0*tE^(5/2)*1.0/nbK





    # println((sum))


    



    return sum 
end 




function _dFdE(E::Float64, Lz::Float64, nbK::Int64=nbK_default)

    return _dFdE_eps(-1.0,E,Lz,nbK) + _dFdE_eps(1.0,E,Lz,nbK)

end



# For now, use the integral definition
# Batsleer & Dejonghe (1993)
# https://ui.adsabs.harvard.edu/abs/1993A%26A...271..104B/abstract
function _dFdLz_eps(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=nbK_default)


    # dF/dtLz
    tE = _tE(E)
    tLz = _tLz(Lz)
    f0 = M/(G*M*(a+c))^(3/2)*(tc^2)/(2.0^(3/2)*pi^3*ta)
    z = 2.0*h*tE*tLz^2

   
    dzdtLz = 2*h*tE*2*tLz

    # Derivative within the integral 
    sum2 = 0.0
    for k=1:nbK
        t = 1.0/nbK*(k-0.5)
        xe = x_eps(tE,eps,z,t)

        dxedtLz =  2*ta*tE*t*sqrt(1-t^2)*(-eps*t*0.5*dzdtLz/sqrt(z))/(1+eps*sqrt(z)*t)^2



        # num = (1 - t^2) * [ (3 + 4 xe - xe^2)(1 - xe)(1-t^2) + 12 t^2 ]
        # den = (1 - 2 ta tE t sqrt(1-t^2) + eps sqrt(z) t)^5

        num = (1.0-t^2)*((3.0+4.0*xe-xe^2)*(1.0-xe)*(1.0-t^2)+12.0*t^2)
        den = (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^5

        dnum = (1.0-t^2)*((4*dxedtLz-2*xe*dxedtLz)*(1-xe)*(1.0-t^2) + (3.0+4.0*xe-xe^2)*(-dxedtLz)*(1.0-t^2))
        dden = 5.0 * ( eps*t*(0.5*dzdtLz/sqrt(z)) ) * (1.0-2.0*ta*tE*t*sqrt(1.0-t^2)+eps*sqrt(z)*t)^4


        sum2 += (dnum*den-num*dden)/den^2 
    end 

    sum2 *= f0*tE^(5/2)*1.0/nbK


    # dF/dLz = dF/dtLz dtLz/dLz
    # dtLz/dELz = 1/sqrt((a+c)*G*M)


    sum2 *= 1/sqrt((a+c)*G*M)

    return sum2
end 

function _dFdLz(E::Float64, Lz::Float64, nbK::Int64=nbK_default)

    return _dFdLz_eps(-1.0,E,Lz,nbK) + _dFdLz_eps(1.0,E,Lz,nbK)

end


















function test_dFdE(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    E_p = E + eps 
    E_m = E - eps 

    dFdE_num = (F(E_p,Lz,nbK) - F(E_m,Lz,nbK))/(2*eps)
    dFdE_th = _dFdE(E,Lz,nbK)

    println("Num = ",dFdE_num)
    println("Th  = ",dFdE_th)

    println("Num/Th =", dFdE_num/dFdE_th)


end


function test_dFdLz(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    Lz_p = Lz + eps 
    Lz_m = Lz - eps 

    dFdLz_num = (F(E,Lz_p,nbK) - F(E,Lz_m,nbK))/(2*eps)
    dFdLz_th = _dFdLz(E,Lz,nbK)

    println("Num = ",dFdLz_num)
    println("Th  = ",dFdLz_th)

    println("Num/Th =", dFdLz_num/dFdLz_th)



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




function test_dFdJ(Ju, Jv, Lz, dJ=0.001, nbK::Int64=100, nbu::Int64=100, nbv::Int64=100, eps::Float64=0.01)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    dJudE, dJudI3, dJudLz = dJu(E,Lz,I3,nbu,eps)
    dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3,nbv)

   detJ = dJudE*dJvdI3-dJudI3*dJvdE

   Omegau = dJvdI3/detJ 
   Omegav = -dJudI3/detJ 
   Omegaz = (dJudI3*dJvdLz-dJudLz*dJvdI3)/detJ 


    dFdE = _dFdE(E,Lz,nbK)
    dFdLz = _dFdLz(E,Lz,nbK)

    dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
    dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
    dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

    println("Theory: ",( dFdJu, dFdJv, dFdLz))

    # Num 

    # Ju
    E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju+dJ,Lz,Jv)
    E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju-dJ,Lz,Jv)

    Fp = F(E_p,Lz_p,nbK)
    Fm = F(E_m,Lz_m,nbK)

    dFdJu = (Fp-Fm)/(2.0*dJ)

     # Jv
     E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv+dJ)
     E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv-dJ)
 
     Fp = F(E_p,Lz_p,nbK)
     Fm = F(E_m,Lz_m,nbK)
 
     dFdJv = (Fp-Fm)/(2.0*dJ)
 
     # Lz
     E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz+dJ,Jv)
     E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz-dJ,Jv)
 
     Fp = F(E_p,Lz_p,nbK)
     Fm = F(E_m,Lz_m,nbK)
 
     dFdLz = (Fp-Fm)/(2.0*dJ)
 
     println("Num   : ",( dFdJu, dFdJv, dFdLz))

end



function test(eps::Float64, E::Float64, Lz::Float64, nbK::Int64=100)

    dE = 0.000001

    z = 2.0*h*E*Lz^2

    sumD = 0.0

    for k=1:nbK 

        t = 1.0/nbK*(k-0.5)

        z = 2*h*E*Lz^2
        xe = 2*a*E*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t)

        dzdE = 2*h*Lz^2
        dxedE = 2*a*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t) - 2*a*E*t*sqrt(1-t^2)*(eps*0.5*dzdE/sqrt(z)*t)/(1+eps*sqrt(z)*t)^2


        num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        den = (1-2*a*E*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

        dnumdE = (4*dxedE-2*dxedE*xe)*(1-xe)*(1-t^2)+  (3+4*xe-xe^2)*(-dxedE)*(1-t^2) 
        ddendE = 5*(1-2*a*E*t*sqrt(1-t^2) + eps*sqrt(z)*t)^4 * (-2*a*t*sqrt(1-t^2) + eps*0.5*dzdE/sqrt(z)*t)

        sumD += (1-t^2)*(dnumdE*den-num*ddendE)/den^2

    end

    sumD *= 1.0/nbK

    println(sumD)


    # num 

    sump = 0.0

    for k=1:nbK 

        t = 1.0/nbK*(k-0.5)

        z = 2*h*(E+dE)*Lz^2
        xe = 2*a*(E+dE)*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t)

        num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        den = (1-2*a*(E+dE)*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

 
        sump += (1-t^2)*num/den 

    end

    sump *= 1.0/nbK

    summ = 0.0

    for k=1:nbK 

        t = 1.0/nbK*(k-0.5)

        z = 2*h*(E-dE)*Lz^2
        xe = 2*a*(E-dE)*t*sqrt(1-t^2)/(1+eps*sqrt(z)*t)

        num = (3+4*xe-xe^2)*(1-xe)*(1-t^2)+ 12*t^2
        den = (1-2*a*(E-dE)*t*sqrt(1-t^2) + eps*sqrt(z)*t)^5

 
        summ += (1-t^2)*num/den 

    end

    summ *= 1.0/nbK

    println((sump-summ)/(2.0*dE))
end


function rho_from_F_test(u::Float64, v::Float64, nbE::Int64=100, nbv::Int64=100, nbK::Int64=100)

    lambda, nu = lambda_nu_from_u_v(u,v)
    R, z = R_z_from_lambda_nu(lambda,nu)
    psiRz = psi(lambda,nu)
    rhoRz = rho(lambda,nu)

    # Exact expression
    println("rho theory = ",rhoRz)

    # Exact with R,z 
    num = (a^2+c^2)*R^2 + 2*a^2*z^2 + 2*a^2*c^2+a^4+3*a^2*sqrt(a^2*c^2+c^2*R^2+a^2*z^2)
    den = (a^2*c^2+c^2*R^2+a^2*z^2)^(3/2)*(R^2+z^2+a^2+c^2+2*sqrt(a^2*c^2+c^2*R^2+a^2*z^2))^(3/2)
    println("rho Rz the = ",M*c^2/(4*pi)*num/den)

    # Integration over F(E,Lz)
    sum = 0.0

    for iE=1:nbE 
        E = psiRz + (0-psiRz)/nbE*(iE-0.5)

        for iv=1:nbv 
            vphi = sqrt(2.0*(E-psiRz))/nbv*(iv-0.5)
            Lz = R*vphi 

            Ftot = F(E,Lz,nbK)

            sum += Ftot*sqrt(2.0*(E-psiRz))/nbv*(0-psiRz)/nbE
        end
    end

    sum *=4.0*pi 

    println("rho integr = ",sum)
end