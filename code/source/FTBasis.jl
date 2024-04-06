# include("RadialBasisExact.jl")
# # include("RadialBasis.jl)
# include("SpheroidalHarmonics.jl")
# include("Args.jl")
# include("Constants.jl")

# using LinearAlgebra


# Use exact formalism here
# https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
# https://github.com/JuliaStellarDynamics/LinearResponse.jl/blob/main/src/WMat.jl
# try to compute all resonance at once ?
# use matrix product property on alphak, betak 
# https://discourse.julialang.org/t/why-putting-partial-array-in-the-function-argument-and-update-it-actually-does-not-change-the-array/66754
function tab_Wkp!(p::Int64, tab_Wkp::Matrix{ComplexF64}, m::Int64, elem::OrbitalElements, tab_WuWv::Matrix{ComplexF64}, tab_ang_k1::Vector{ComplexF64}, tab_ang_k2::Vector{ComplexF64},  matrix_grad::Matrix{Float64}, nbk::Int64, nbt::Int64=nbt_default)

    # k1,k2 in [-kmax,kmax]
    # k3 = m
    # everthing in elem


    # Compute alpha_k, beta_k on the fly 
    # Backward integration

    # tab_Wkp[k]
    # tab_Wkp[k,i] with i=1,2,3,4 for u1, u2, v1, v2

    # tab_alphak[k] : alphak
    # tab_betak[k] : betak

    Jw = abs(matrix_grad[1,1]*matrix_grad[2,2]-matrix_grad[2,1]*matrix_grad[1,2])

    for k=1:nbk
        for i=1:4
            tab_WuWv[k,i] = 0.0 + 0.0im 

        end


    end

    # alphak = k1 alpha1 + k2 alpha2 + m alpha3
    alpha_1, alpha_2, alpha_3 = pi, 0.0, 0.0

  

    # Backward integration 
    # RK4

    du = pi/nbt 


    iu = 1 # index to access elem.tab_ylm
    sinhu = elem.tab_sinhu[iu]
    tpu = elem.tab_tpu[iu]
    Jwu = Delta^4/Jw * sinhu^2

    Flmn = elem.tab_Flmn[iu]

    

    for i=1:nbt 

        # Step 1


        # alphak = k1 alpha1 

        # Start at k1=-kmax, k2=-kmax 
        # exp (-i alpha_k) = exp(- i k_1 alpha_1) exp(- i k_2 alpha_2) exp(- i m alpha_3)
        # k_1, k_2 from -kmax to kmax 
        # Start at k_1, k_2 = kmax with exp( i k_max alpha_1) and exp( i k_max alpha_2)
        # Then multiply each term by exp(- i alpha_1) and exp(- i alpha_2)

        exp_minus_i_alphak_3 = exp(-1.0im*m*alpha_3)

        exp_minus_i_alpha1 = exp(-1.0im*alpha_1)
        exp_minus_i_alpha2 = exp(-1.0im*alpha_2)
        
        # pref_Wu2_1 = du/6.0 * 1.0/pi * Flmn/tpu *exp_minus_i_alphak_3
        pref_Wu2_1 = 1.0/6.0 * Flmn/tpu *exp_minus_i_alphak_3
        pref_Wu1_1 = Jwu * pref_Wu2_1


        # compute k and -k at the same time using properties of exp

        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_alpha1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_alpha2
        end


        for k1=0:kmax 

            ang1 = tab_ang_k1[k1+1]

            # prefu1_k1 = pref_Wu1_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3
            # prefu2_k1 = pref_Wu2_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3

            prefu1_k1 = pref_Wu1_1*ang1 
            prefu2_k1 = pref_Wu2_1*ang1 

            # prefu1_minus_k1 = pref_Wu1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3
            # prefu2_minus_k1 = pref_Wu2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3

            prefu1_minus_k1 = pref_Wu1_1*conj(ang1)
            prefu2_minus_k1 = pref_Wu2_1*conj(ang1)


            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,1] += prefu1_k1 * ang2
                tab_WuWv[k,2] += prefu2_k1 * ang2   

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,1] += prefu1_k1 * conj(ang2)
                    tab_WuWv[k,2] += prefu2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,1] += prefu1_minus_k1 * ang2
                    tab_WuWv[k,2] += prefu2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,1] += prefu1_minus_k1 * conj(ang2)
                        tab_WuWv[k,2] += prefu2_minus_k1 * conj(ang2)
                    end
                end
            end
        end


        # Update alpha's

        # dalpha_E_1 = du*(matrix_freq[1,1]*elem.tab_alpha_E[iu] +  matrix_freq[2,1]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,1]*elem.tab_alpha_Lz[iu] )
        # dalpha_I3_1 = du*(matrix_freq[1,2]*elem.tab_alpha_E[iu] +  matrix_freq[2,2]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,2]*elem.tab_alpha_Lz[iu] )
        # dalpha_Lz_1 = du*(matrix_freq[1,3]*elem.tab_alpha_E[iu] +  matrix_freq[2,3]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,3]*elem.tab_alpha_Lz[iu] )


        dalpha_E_1 = du*(elem.tab_alpha_E[iu] )
        dalpha_I3_1 = du*( elem.tab_alpha_I3[iu]    )
        dalpha_Lz_1 = du*(elem.tab_alpha_Lz[iu] )




         # Step 2

         iu += 1 # index to access elem.tab_ylm
         sinhu = elem.tab_sinhu[iu]
         tpu = elem.tab_tpu[iu]
         Jwu = Delta^4/Jw * sinhu^2
 
         Flmn = elem.tab_Flmn[iu]

         exp_minus_i_alphak_3 = exp(-1.0im*m*(alpha_3-0.5*dalpha_Lz_1))

         exp_minus_i_alpha1 = exp(-1.0im*(alpha_1-0.5*dalpha_E_1))
         exp_minus_i_alpha2 = exp(-1.0im*(alpha_2-0.5*dalpha_I3_1))
         
        #  pref_Wu2_2 = du/3.0 * 1.0/pi * Flmn/tpu*exp_minus_i_alphak_3
         pref_Wu2_2 = 1.0/3.0   * Flmn/tpu*exp_minus_i_alphak_3
         pref_Wu1_2 = Jwu * pref_Wu2_2



        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_alpha1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_alpha2
        end


        for k1=0:kmax 


            ang1 = tab_ang_k1[k1+1]

            # prefu1_k1 = pref_Wu1_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3
            # prefu2_k1 = pref_Wu2_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3

            prefu1_k1 = pref_Wu1_2*ang1 
            prefu2_k1 = pref_Wu2_2*ang1 

            # prefu1_minus_k1 = pref_Wu1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3
            # prefu2_minus_k1 = pref_Wu2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3

            prefu1_minus_k1 = pref_Wu1_2*conj(ang1)
            prefu2_minus_k1 = pref_Wu2_2*conj(ang1)

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]

                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,1] += prefu1_k1 * ang2
                tab_WuWv[k,2] += prefu2_k1 * ang2

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,1] += prefu1_k1 * conj(ang2)
                    tab_WuWv[k,2] += prefu2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,1] += prefu1_minus_k1 * ang2
                    tab_WuWv[k,2] += prefu2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,1] += prefu1_minus_k1 * conj(ang2)
                        tab_WuWv[k,2] += prefu2_minus_k1 * conj(ang2)
                    end
                end
            end
        end



        # Update alphak 
  

        # dalpha_E_2 = du*(matrix_freq[1,1]*elem.tab_alpha_E[iu] +  matrix_freq[2,1]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,1]*elem.tab_alpha_Lz[iu] )
        # dalpha_I3_2 = du*(matrix_freq[1,2]*elem.tab_alpha_E[iu] +  matrix_freq[2,2]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,2]*elem.tab_alpha_Lz[iu] )
        # dalpha_Lz_2 = du*(matrix_freq[1,3]*elem.tab_alpha_E[iu] +  matrix_freq[2,3]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,3]*elem.tab_alpha_Lz[iu] )


        dalpha_E_2 = du*(elem.tab_alpha_E[iu] )
        dalpha_I3_2 = du*( elem.tab_alpha_I3[iu]    )
        dalpha_Lz_2 = du*(elem.tab_alpha_Lz[iu] )

        # Step 3


        exp_minus_i_alphak_3 = exp(-1.0im*m*(alpha_3-0.5*dalpha_Lz_2))

        exp_minus_i_alpha1 = exp(-1.0im*(alpha_1-0.5*dalpha_E_2))
        exp_minus_i_alpha2 = exp(-1.0im*(alpha_2-0.5*dalpha_I3_2))
        
        # pref_Wu2_3 = du/3.0 * 1.0/pi * Flmn/tpu*exp_minus_i_alphak_3
        pref_Wu2_3 = 1.0/3.0  * Flmn/tpu*exp_minus_i_alphak_3
        pref_Wu1_3 = Jwu * pref_Wu2_3



        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_alpha1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_alpha2
        end


        for k1=0:kmax 


            ang1 = tab_ang_k1[k1+1]

            # prefu1_k1 = pref_Wu1_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3
            # prefu2_k1 = pref_Wu2_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3

            prefu1_k1 = pref_Wu1_3*ang1 
            prefu2_k1 = pref_Wu2_3*ang1 

            # prefu1_minus_k1 = pref_Wu1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3
            # prefu2_minus_k1 = pref_Wu2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3

            prefu1_minus_k1 = pref_Wu1_3*conj(ang1)
            prefu2_minus_k1 = pref_Wu2_3*conj(ang1)

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,1] += prefu1_k1 * ang2
                tab_WuWv[k,2] += prefu2_k1 * ang2

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,1] += prefu1_k1 * conj(ang2)
                    tab_WuWv[k,2] += prefu2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,1] += prefu1_minus_k1 * ang2
                    tab_WuWv[k,2] += prefu2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,1] += prefu1_minus_k1 * conj(ang2)
                        tab_WuWv[k,2] += prefu2_minus_k1 * conj(ang2)
                    end
                end
            end
        end


        # Update alphak 
  

        dalpha_E_3 = dalpha_E_2
        dalpha_I3_3 = dalpha_I3_2
        dalpha_Lz_3 =  dalpha_Lz_2

        # Step 4

        iu += 1 # index to access elem.tab_ylm
        sinhu = elem.tab_sinhu[iu]
        tpu = elem.tab_tpu[iu]
        Jwu = Delta^4/Jw * sinhu^2

        Flmn = elem.tab_Flmn[iu]

        exp_minus_i_alphak_3 = exp(-1.0im*m*(alpha_3-dalpha_Lz_3))

        exp_minus_i_alpha1 = exp(-1.0im*(alpha_1-dalpha_E_3))
        exp_minus_i_alpha2 = exp(-1.0im*(alpha_2-dalpha_I3_3))
        
        # pref_Wu2_4 = du/6.0 * 1.0/pi * Flmn/tpu*exp_minus_i_alphak_3
        pref_Wu2_4 = 1.0/6.0   * Flmn/tpu*exp_minus_i_alphak_3
        pref_Wu1_4 = Jwu * pref_Wu2_4

   

        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_alpha1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_alpha2
        end


        for k1=0:kmax 


            ang1 = tab_ang_k1[k1+1]

            # prefu1_k1 = pref_Wu1_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3
            # prefu2_k1 = pref_Wu2_1*tab_ang_k1[k1+1]#*exp_minus_i_alphak_3

            prefu1_k1 = pref_Wu1_4*ang1 
            prefu2_k1 = pref_Wu2_4*ang1 

            # prefu1_minus_k1 = pref_Wu1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3
            # prefu2_minus_k1 = pref_Wu2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_alphak_3

            prefu1_minus_k1 = pref_Wu1_4*conj(ang1)
            prefu2_minus_k1 = pref_Wu2_4*conj(ang1)

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,1] += prefu1_k1 * ang2
                tab_WuWv[k,2] += prefu2_k1 * ang2  

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,1] += prefu1_k1 * conj(ang2)
                    tab_WuWv[k,2] += prefu2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,1] += prefu1_minus_k1 * ang2
                    tab_WuWv[k,2] += prefu2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,1] += prefu1_minus_k1 * conj(ang2)
                        tab_WuWv[k,2] += prefu2_minus_k1 * conj(ang2)
                    end
                end
            end
        end

        # Update alphak 
        
       # dalpha_E_4 = du*(matrix_freq[1,1]*elem.tab_alpha_E[iu] +  matrix_freq[2,1]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,1]*elem.tab_alpha_Lz[iu] )
        # dalpha_I3_4 = du*(matrix_freq[1,2]*elem.tab_alpha_E[iu] +  matrix_freq[2,2]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,2]*elem.tab_alpha_Lz[iu] )
        # dalpha_Lz_4 = du*(matrix_freq[1,3]*elem.tab_alpha_E[iu] +  matrix_freq[2,3]*elem.tab_alpha_I3[iu]  +  matrix_freq[3,3]*elem.tab_alpha_Lz[iu] )


        dalpha_E_4 = du*(elem.tab_alpha_E[iu] )
        dalpha_I3_4 = du*( elem.tab_alpha_I3[iu]    )
        dalpha_Lz_4 = du*(elem.tab_alpha_Lz[iu] )
        
        
        # Update alphak for next iteration
        alpha_1 -= (dalpha_E_1 + 2.0*dalpha_E_2 + 2.0*dalpha_E_3 + dalpha_E_4)/6.0
        alpha_2 -= (dalpha_I3_1 + 2.0*dalpha_I3_2 + 2.0*dalpha_I3_3 + dalpha_I3_4)/6.0
        alpha_3 -= (dalpha_Lz_1 + 2.0*dalpha_Lz_2 + 2.0*dalpha_Lz_3 + dalpha_Lz_4)/6.0

        # alpha_1 -= (dalpha_E_1 + 4.0*dalpha_E_2+ dalpha_E_4)/6.0
        # alpha_2 -= (dalpha_I3_1 + 4.0*dalpha_I3_2 + dalpha_I3_4)/6.0
        # alpha_3 -= (dalpha_Lz_1 + 4.0*dalpha_Lz_2 + dalpha_Lz_4)/6.0


    end



    # BETA

        
    # betak = k1 beta1 + k2 beta2 + m beta3
    beta_1, beta_2, beta_3 = 0.0, pi, 0.0

    dv = pi/nbt 

    iv = 1 # index to access elem.tab_ylm
    sinv = elem.tab_sinv[iv]
    ylm = elem.tab_ylm[iv]
    tpv = elem.tab_tpv[iv]
    Jwv = Delta^4/Jw * sinv^2
    


     for i=1:nbt 

        # Step 1

        exp_minus_i_betak_3 = exp(-1.0im*m*beta_3)

        exp_minus_i_beta1 = exp(-1.0im*beta_1)
        exp_minus_i_beta2 = exp(-1.0im*beta_2)
        
        # pref_Wv1_1 = dv/6.0 * 1/pi * ylm/tpv*exp_minus_i_betak_3
        pref_Wv1_1 = 1.0/6.0  * ylm/tpv*exp_minus_i_betak_3
        pref_Wv2_1 = Jwv * pref_Wv1_1

     

        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_beta1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_beta2
        end

        for k1=0:kmax 

            ang1 = tab_ang_k1[k1+1]

            # prefv1_k1 = pref_Wv1_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3
            # prefv2_k1 = pref_Wv2_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3

            prefv1_k1 = pref_Wv1_1*ang1#*exp_minus_i_betak_3
            prefv2_k1 = pref_Wv2_1*ang1#*exp_minus_i_betak_3

            # prefv1_minus_k1 = pref_Wv1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3
            # prefv2_minus_k1 = pref_Wv2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3

            prefv1_minus_k1 = pref_Wv1_1*conj(ang1)#*exp_minus_i_betak_3
            prefv2_minus_k1 = pref_Wv2_1*conj(ang1)#*exp_minus_i_betak_3

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,3] += prefv1_k1 * ang2
                tab_WuWv[k,4] += prefv2_k1 * ang2  

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,3] += prefv1_k1 * conj(ang2)
                    tab_WuWv[k,4] += prefv2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,3] += prefv1_minus_k1 * ang2
                    tab_WuWv[k,4] += prefv2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,3] += prefv1_minus_k1 * conj(ang2)
                        tab_WuWv[k,4] += prefv2_minus_k1 * conj(ang2)
                    end
                end
            end
        end

        # Update betak 

        # dbeta_E_1 = du*(matrix_freq[1,1]*elem.tab_beta_E[iv] +  matrix_freq[2,1]*elem.tab_beta_I3[iv]  +  matrix_freq[3,1]*elem.tab_beta_Lz[iv] )
        # dbeta_I3_1 = du*(matrix_freq[1,2]*elem.tab_beta_E[iv] +  matrix_freq[2,2]*elem.tab_beta_I3[iv]  +  matrix_freq[3,2]*elem.tab_beta_Lz[iv] )
        # dbeta_Lz_1 = du*(matrix_freq[1,3]*elem.tab_beta_E[iv] +  matrix_freq[2,3]*elem.tab_beta_I3[iv]  +  matrix_freq[3,3]*elem.tab_beta_Lz[iv] )

        dbeta_E_1 = du*(elem.tab_beta_E[iv]  )
        dbeta_I3_1 = du*(elem.tab_beta_I3[iv]  )
        dbeta_Lz_1 = du*(elem.tab_beta_Lz[iv] )


         # Step 2

         iv += 1 # update
         sinv = elem.tab_sinv[iv]
         tpv = elem.tab_tpv[iv]
         Jwv = Delta^4/Jw * sinv^2
         ylm = elem.tab_ylm[iv]

        exp_minus_i_betak_3 = exp(-1.0im*m*(beta_3-0.5*dbeta_Lz_1))

        exp_minus_i_beta1 = exp(-1.0im*(beta_1-0.5*dbeta_E_1))
        exp_minus_i_beta2 = exp(-1.0im*(beta_2-0.5*dbeta_I3_1))
        
        # pref_Wv1_2 = dv/3.0 * 1/pi * ylm/tpv*exp_minus_i_betak_3
        pref_Wv1_2 = 1.0/3.0 * ylm/tpv*exp_minus_i_betak_3
        pref_Wv2_2 = Jwv * pref_Wv1_2

    
        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_beta1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_beta2
        end

        for k1=0:kmax 

            ang1 = tab_ang_k1[k1+1]

            # prefv1_k1 = pref_Wv1_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3
            # prefv2_k1 = pref_Wv2_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3

            prefv1_k1 = pref_Wv1_2*ang1#*exp_minus_i_betak_3
            prefv2_k1 = pref_Wv2_2*ang1#*exp_minus_i_betak_3

            # prefv1_minus_k1 = pref_Wv1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3
            # prefv2_minus_k1 = pref_Wv2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3

            prefv1_minus_k1 = pref_Wv1_2*conj(ang1)#*exp_minus_i_betak_3
            prefv2_minus_k1 = pref_Wv2_2*conj(ang1)#*exp_minus_i_betak_3

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,3] += prefv1_k1 * ang2
                tab_WuWv[k,4] += prefv2_k1 * ang2 

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,3] += prefv1_k1 * conj(ang2)
                    tab_WuWv[k,4] += prefv2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,3] += prefv1_minus_k1 * ang2
                    tab_WuWv[k,4] += prefv2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,3] += prefv1_minus_k1 * conj(ang2)
                        tab_WuWv[k,4] += prefv2_minus_k1 * conj(ang2)
                    end
                end
            end
        end

        # Update betak 

        # dbeta_E_2 = du*(matrix_freq[1,1]*elem.tab_beta_E[iv] +  matrix_freq[2,1]*elem.tab_beta_I3[iv]  +  matrix_freq[3,1]*elem.tab_beta_Lz[iv] )
        # dbeta_I3_2 = du*(matrix_freq[1,2]*elem.tab_beta_E[iv] +  matrix_freq[2,2]*elem.tab_beta_I3[iv]  +  matrix_freq[3,2]*elem.tab_beta_Lz[iv] )
        # dbeta_Lz_2 = du*(matrix_freq[1,3]*elem.tab_beta_E[iv] +  matrix_freq[2,3]*elem.tab_beta_I3[iv]  +  matrix_freq[3,3]*elem.tab_beta_Lz[iv] )

        dbeta_E_2 = du*(elem.tab_beta_E[iv]  )
        dbeta_I3_2 = du*(elem.tab_beta_I3[iv]  )
        dbeta_Lz_2 = du*(elem.tab_beta_Lz[iv] )


         # Step 3

        exp_minus_i_betak_3 = exp(-1.0im*m*(beta_3-0.5*dbeta_Lz_2))

        exp_minus_i_beta1 = exp(-1.0im*(beta_1-0.5*dbeta_E_2))
        exp_minus_i_beta2 = exp(-1.0im*(beta_2-0.5*dbeta_I3_2))
        
        # pref_Wv1_3 = dv/3.0 * 1/pi * ylm/tpv*exp_minus_i_betak_3
        pref_Wv1_3 = 1.0/3.0  * ylm/tpv*exp_minus_i_betak_3
        pref_Wv2_3 = Jwv * pref_Wv1_3



        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_beta1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_beta2
        end

        for k1=0:kmax 

            ang1 = tab_ang_k1[k1+1]

            # prefv1_k1 = pref_Wv1_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3
            # prefv2_k1 = pref_Wv2_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3

            prefv1_k1 = pref_Wv1_3*ang1#*exp_minus_i_betak_3
            prefv2_k1 = pref_Wv2_3*ang1#*exp_minus_i_betak_3

            # prefv1_minus_k1 = pref_Wv1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3
            # prefv2_minus_k1 = pref_Wv2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3

            prefv1_minus_k1 = pref_Wv1_3*conj(ang1)#*exp_minus_i_betak_3
            prefv2_minus_k1 = pref_Wv2_3*conj(ang1)#*exp_minus_i_betak_3

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,3] += prefv1_k1 * ang2
                tab_WuWv[k,4] += prefv2_k1 * ang2      

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,3] += prefv1_k1 * conj(ang2)
                    tab_WuWv[k,4] += prefv2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,3] += prefv1_minus_k1 * ang2
                    tab_WuWv[k,4] += prefv2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,3] += prefv1_minus_k1 * conj(ang2)
                        tab_WuWv[k,4] += prefv2_minus_k1 * conj(ang2)
                    end
                end
            end
        end

        # Update betak 
  

        dbeta_E_3 = dbeta_E_2
        dbeta_I3_3 = dbeta_I3_2
        dbeta_Lz_3 =  dbeta_Lz_2


         # Step 4

         iv += 1 # update
         sinv = elem.tab_sinv[iv]
         tpv = elem.tab_tpv[iv]
         Jwv = Delta^4/Jw * sinv^2
         ylm = elem.tab_ylm[iv]

        exp_minus_i_betak_3 = exp(-1.0im*m*(beta_3-dbeta_Lz_3))

        exp_minus_i_beta1 = exp(-1.0im*(beta_1-dbeta_E_3))
        exp_minus_i_beta2 = exp(-1.0im*(beta_2-dbeta_I3_3))
        
        # pref_Wv1_4 = dv/6.0 * 1/pi * ylm/tpv*exp_minus_i_betak_3
        pref_Wv1_4 = 1.0/6.0 * ylm/tpv*exp_minus_i_betak_3
        pref_Wv2_4 = Jwv * pref_Wv1_4

    
      
        tab_ang_k1[1] = 1.0 + 0.0im # Initialize at k1=0 
        tab_ang_k2[1] = 1.0 + 0.0im # Initialize at k2=0

        for k=1:kmax
            tab_ang_k1[k+1] =  tab_ang_k1[k]*exp_minus_i_beta1
            tab_ang_k2[k+1] =  tab_ang_k2[k]*exp_minus_i_beta2
        end

        for k1=0:kmax 

            ang1 = tab_ang_k1[k1+1]

            # prefv1_k1 = pref_Wv1_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3
            # prefv2_k1 = pref_Wv2_1*tab_ang_k1[k1+1]#*exp_minus_i_betak_3

            prefv1_k1 = pref_Wv1_4*ang1#*exp_minus_i_betak_3
            prefv2_k1 = pref_Wv2_4*ang1#*exp_minus_i_betak_3

            # prefv1_minus_k1 = pref_Wv1_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3
            # prefv2_minus_k1 = pref_Wv2_1*conj(tab_ang_k1[k1+1])#*exp_minus_i_betak_3

            prefv1_minus_k1 = pref_Wv1_4*conj(ang1)#*exp_minus_i_betak_3
            prefv2_minus_k1 = pref_Wv2_4*conj(ang1)#*exp_minus_i_betak_3

            for k2=0:kmax 

                ang2 = tab_ang_k2[k2+1]


                # +k1 +k2
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                tab_WuWv[k,3] += prefv1_k1 * ang2
                tab_WuWv[k,4] += prefv2_k1 * ang2   

                if (k2 != 0)
                    # +k1 -k2
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    tab_WuWv[k,3] += prefv1_k1 * conj(ang2)
                    tab_WuWv[k,4] += prefv2_k1 * conj(ang2)
                end
                 
                if (k1 != 0)
                    # -k1 +k2
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    tab_WuWv[k,3] += prefv1_minus_k1 * ang2
                    tab_WuWv[k,4] += prefv2_minus_k1 * ang2

                    if (k2 != 0)

                        # -k1 -k2
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                        tab_WuWv[k,3] += prefv1_minus_k1 * conj(ang2)
                        tab_WuWv[k,4] += prefv2_minus_k1 * conj(ang2)
                    end
                end
            end
        end

        # Update betak 

         # dbeta_E_4 = du*(matrix_freq[1,1]*elem.tab_beta_E[iv] +  matrix_freq[2,1]*elem.tab_beta_I3[iv]  +  matrix_freq[3,1]*elem.tab_beta_Lz[iv] )
        # dbeta_I3_4 = du*(matrix_freq[1,2]*elem.tab_beta_E[iv] +  matrix_freq[2,2]*elem.tab_beta_I3[iv]  +  matrix_freq[3,2]*elem.tab_beta_Lz[iv] )
        # dbeta_Lz_4 = du*(matrix_freq[1,3]*elem.tab_beta_E[iv] +  matrix_freq[2,3]*elem.tab_beta_I3[iv]  +  matrix_freq[3,3]*elem.tab_beta_Lz[iv] )

        dbeta_E_4 = du*(elem.tab_beta_E[iv]  )
        dbeta_I3_4 = du*(elem.tab_beta_I3[iv]  )
        dbeta_Lz_4 = du*(elem.tab_beta_Lz[iv] )


        # Update betak for next iteration
        beta_1 -= (dbeta_E_1 + 2.0*dbeta_E_2 + 2.0*dbeta_E_3 + dbeta_E_4)/6.0
        beta_2 -= (dbeta_I3_1 + 2.0*dbeta_I3_2 + 2.0*dbeta_I3_3 + dbeta_I3_4)/6.0
        beta_3 -= (dbeta_Lz_1 + 2.0*dbeta_Lz_2 + 2.0*dbeta_Lz_3 + dbeta_Lz_4)/6.0

        # beta_1 -= (dbeta_E_1 + 4.0*dbeta_E_2 + dbeta_E_4)/6.0
        # beta_2 -= (dbeta_I3_1 + 4.0*dbeta_I3_2 + dbeta_I3_4)/6.0
        # beta_3 -= (dbeta_Lz_1 + 4.0*dbeta_Lz_2  + dbeta_Lz_4)/6.0



    end


    for k=1:nbk 

        # du/pi = dv/pi = 1/nbt
        tab_Wkp[p,k] = (tab_WuWv[k,1]*tab_WuWv[k,3] + tab_WuWv[k,2]*tab_WuWv[k,4])/(nbt*nbt)

    end

    # println("===")
    # display(tab_Wkp[p,:])

end
 



# Optimize for frequency sampling
function ResponseMatrix_m_sampling(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

    @assert (abs(m) <= mmax) "m must be less that m_max"

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbgrid = nbJ^3
    nbk = (2*kmax+1)^2
    nbomega = nbRe * nbIm 

    dJu = Jumax/nbJ
    dJv = Jvmax/nbJ
    dLz = 2.0*Lzmax/nbJ

    tab_omega = zeros(Float64, nbomega, 2) # (re, im)

    for iomega=1:nbomega 

        # iomega - 1 = im-1 + nbIm*(ire-1)

        ire = floor(Int64,(iomega-1)/nbIm) + 1
        iim = iomega - nbIm*(ire-1)

        re_omega = re_omega_min + (re_omega_max-re_omega_min)/(nbRe-1)*(ire-1)
        im_omega = im_omega_min + (im_omega_max-im_omega_min)/(nbIm-1)*(iim-1)

        tab_omega[iomega,1], tab_omega[iomega,2] = re_omega, im_omega 
    end


    elem = [OrbitalElements_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

    # Compute table of Mpq

    # tab_Mpq_par_re = [Threads.Atomic{Float64}(0.0) for index=1:nbp*nbp, iomega=1:nbomega] 
    # tab_Mpq_par_im = [Threads.Atomic{Float64}(0.0) for index=1:nbp*nbp, iomega=1:nbomega] 

    # tab_Mpq_par_re_par = [0.0 for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()] 
    # tab_Mpq_par_im_par = [0.0 for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()] 

    tab_Mpq_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()] 

    tabWkp_temp = [zeros(ComplexF64, nbp, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(ComplexF64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_ang_k1 = [zeros(ComplexF64, kmax+1) for it=1:Threads.nthreads()] 
    tab_ang_k2 = [zeros(ComplexF64, kmax+1) for it=1:Threads.nthreads()] 
    
    # Fill the tables of Wkp

    Threads.@threads for iJ=1:nbgrid 

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ 
        # nbv = nbJ 

        iu = floor(Int64,(iJ - 1)/nbJ^2) + 1
        leftover = iJ - 1 - nbJ^2*(iu-1)

        # leftover = iz-1 + nbz*(iv-1)
        # nbz = nbJ 

        iv = floor(Int64,leftover/nbJ) + 1
        iz = leftover - nbJ*(iv-1) + 1

        # Ju = Jumax/nbJ*(iu-0.5)
        # Jv = Jvmax/nbJ*(iv-0.5)
        # Lz = -Lzmax + 2.0*Lzmax/nbJ*(iz-0.5)

        Ju = dJu*(iu-0.5)
        Jv = dJv*(iv-0.5)
        Lz = -Lzmax + dLz*(iz-0.5)


        E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

        u0, u1 = find_bounds_u(E,Lz,I3)
        v0, v1 = find_bounds_v(E,Lz,I3)

        id_thread = Threads.threadid()

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E,Lz,I3)
        # trans_freq = transpose(freq_matrix[id_thread])
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  



        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        OrbitalElements_update!(elem[id_thread],E,I3,Lz,u0,u1,v0,v1,nbt,freq_matrix[id_thread]) # l does not change: We don't have to update the spheroidal harmonics


        for p=1:nbp 

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l 

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))
            
            

            if (n==n0)
                OrbitalElements_update_ylm!(elem[id_thread],l,m,nbt)
            end

            OrbitalElements_update_Flmn!(elem[id_thread],l,m,n,nbt)
           

            
            tab_Wkp!(p,tabWkp_temp[id_thread],m,elem[id_thread],tabWuWv[id_thread],tab_ang_k1[id_thread],tab_ang_k2[id_thread],grad_matrix[id_thread],nbk,nbt)
           
        end


       for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)
            
            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            kdotdF = k1*dFdJu + k2*dFdJv + m*dFdLz
                        
            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz 

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = kdotdF/(omega-kdotOmega)

                for index=1:nbp*nbp 

                    # index - 1 = (q-1) + (p-1)*nbp
    
                    p = floor(Int64, (index-1)/nbp) + 1
                    q  = index  - (p-1)*nbp
                            
    
                    Wkp = tabWkp_temp[id_thread][p,k]  
                    Wkq = tabWkp_temp[id_thread][q,k]  
    
                    integrand = conj(Wkp) * Wkq * invden

                    # Threads.atomic_add!(tab_Mpq_par_re[index,iomega],real(integrand))
                    # Threads.atomic_add!(tab_Mpq_par_im[index,iomega],imag(integrand))

                    # tab_Mpq_par_re_par[index,iomega,id_thread] += real(integrand)
                    # tab_Mpq_par_im_par[index,iomega,id_thread] += imag(integrand)

                    tab_Mpq_par[index,iomega,id_thread] += integrand
                end

            end

        end

        # println("-----")

        

    end

    tab_Mpq = zeros(ComplexF64, nbp, nbp, nbomega)

    # BE CAREFUL: THERE WAS A MISTAKE HERE BEFORE
    # sqrt(4 pi G)/Delta COMES FROM BOTH Wkp 
    # THEREFORE, WE MUST MULTIPLY BY ITS SQUARE: (4 pi G)/Delta^2
    # pref = (2*pi)^3 * sqrt(4*pi*G)/Delta * Jumax/nbJ * Jvmax/nbJ * 2.0*Lzmax/nbJ
    pref = (2*pi)^3 * 4*pi*G/Delta^2 * Jumax/nbJ * Jvmax/nbJ * 2.0*Lzmax/nbJ

    Threads.@threads for index=1:nbp*nbp 

        # index - 1 = (q-1) + (p-1)*nbp

        p = floor(Int64, (index-1)/nbp) + 1
        q  = index  - (p-1)*nbp

        for iomega=1:nbomega
            

            # tab_Mpq[p,q,iomega] = tab_Mpq_par_re[index,iomega][] + 1im*tab_Mpq_par_im[index,iomega][]

            for it=1:Threads.nthreads()

                # tab_Mpq[p,q,iomega] += tab_Mpq_par_re_par[index,iomega,it] + 1im*tab_Mpq_par_im_par[index,iomega,it]
                tab_Mpq[p,q,iomega] +=tab_Mpq_par[index,iomega,it]

            end
            tab_Mpq[p,q,iomega] *= pref
        end
    end

    return tab_Mpq, tab_omega
end








function det_Dielectric(tabMpq, size)

    id = 1.0* Matrix(I, size, size)

    return det(id - tabMpq)

end





