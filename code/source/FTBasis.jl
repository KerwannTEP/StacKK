
function increment_Wu_half!(m::Int64, l::Int64, tabWuWv::Matrix{Float64}, pref::Float64, Jacu::Float64, alpha_1::Float64, alpha_2::Float64, alpha_3::Float64) #, tab_ang_k1::Array{ComplexF64}, tab_ang_k2::Array{ComplexF64})

    sinalpha1, cosalpha1 = sincos(alpha_1)
    sinalpha2, cosalpha2 = sincos(alpha_2)
    sinmalpha3, cosmalpha3 = sincos(m*alpha_3)


    cosalphak1 = 1.0
    sinalphak1 = 0.0

    for k1=0:kmax 
        cosalphak2 = 1.0
        sinalphak2 = 0.0

        for k2=0:kmax

            if (mod(l+m+k2,2)==0)
            
                # k1 >= 0 ; k2 >= 0
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)

                cosalphak = cosalphak1*cosalphak2*cosmalpha3 - sinalphak1*sinalphak2*cosmalpha3 - sinalphak1*cosalphak2*sinmalpha3 - cosalphak1*sinalphak2*sinmalpha3
                inc = pref * cosalphak
                tabWuWv[k,1] += inc
                tabWuWv[k,2] += inc * Jacu 

                if (k2 != 0)

                    # k1 >= 0 ; k2 <= 0
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)

                    cosalphak = cosalphak1*cosalphak2*cosmalpha3 + sinalphak1*sinalphak2*cosmalpha3 - sinalphak1*cosalphak2*sinmalpha3 + cosalphak1*sinalphak2*sinmalpha3
                    inc = pref * cosalphak
                    tabWuWv[k,1] += inc
                    tabWuWv[k,2] += inc * Jacu 

                end

                if (k1 != 0)
                    # k1 <= 0 ; k2 >= 0
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)

                    cosalphak = cosalphak1*cosalphak2*cosmalpha3 + sinalphak1*sinalphak2*cosmalpha3 + sinalphak1*cosalphak2*sinmalpha3 - cosalphak1*sinalphak2*sinmalpha3
                    inc = pref * cosalphak
                    tabWuWv[k,1] += inc
                    tabWuWv[k,2] += inc * Jacu  

                    if (k2 != 0)

                        # k1 <= 0 ; k2 <= 0
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)
                        # alphak = -k1*alpha_1 - k2*alpha_2 + m*alpha_3 
                        cosalphak = cosalphak1*cosalphak2*cosmalpha3 - sinalphak1*sinalphak2*cosmalpha3 + sinalphak1*cosalphak2*sinmalpha3 + cosalphak1*sinalphak2*sinmalpha3
                        inc = pref * cosalphak
                        tabWuWv[k,1] += inc
                        tabWuWv[k,2] += inc * Jacu 

                    end

                end
            end

            cosalphak2, sinalphak2 = cosalphak2*cosalpha2 - sinalphak2*sinalpha2, sinalphak2*cosalpha2 + cosalphak2*sinalpha2
        end
        cosalphak1, sinalphak1 = cosalphak1*cosalpha1 - sinalphak1*sinalpha1, sinalphak1*cosalpha1 + cosalphak1*sinalpha1
 
    end

end


function increment_Wv_half!(m::Int64, l::Int64, tabWuWv::Matrix{Float64}, pref::Float64, Jacv::Float64, beta_1::Float64, beta_2::Float64, beta_3::Float64)

    sinbeta1, cosbeta1 = sincos(beta_1)
    sinbeta2, cosbeta2 = sincos(beta_2)
    sinmbeta3, cosmbeta3 = sincos(m*beta_3)

    # for k=1:kmax
    #     tab_ang_k1[k+1] =  tab_ang_k1[k]*i_alpha1
    #     tab_ang_k2[k+1] =  tab_ang_k2[k]*i_alpha2
    # end

    cosbetak1 = 1.0
    sinbetak1 = 0.0

    for k1=0:kmax 
        cosbetak2 = 1.0
        sinbetak2 = 0.0
        for k2=0:kmax

            if (mod(l+m+k2,2)==0)

                # k1 >= 0 ; k2 >= 0
                k = 1 + k2 + kmax + (2*kmax+1) * (k1 + kmax)
                # betak = k1*beta_1 + k2*beta_2 + m*beta_3 
                cosbetak = cosbetak1*cosbetak2*cosmbeta3 - sinbetak1*sinbetak2*cosmbeta3 - sinbetak1*cosbetak2*sinmbeta3 - cosbetak1*sinbetak2*sinmbeta3
                inc = pref*cosbetak
                tabWuWv[k,3] += inc * Jacv 
                tabWuWv[k,4] += inc

                if (k2 != 0)

                    # k1 >= 0 ; k2 <= 0
                    k = 1 - k2 + kmax + (2*kmax+1) * (k1 + kmax)
                    # betak = k1*beta_1 - k2*beta_2 + m*beta_3 
                    cosbetak = cosbetak1*cosbetak2*cosmbeta3 + sinbetak1*sinbetak2*cosmbeta3 - sinbetak1*cosbetak2*sinmbeta3 + cosbetak1*sinbetak2*sinmbeta3
                    inc = pref*cosbetak
                    tabWuWv[k,3] += inc * Jacv 
                    tabWuWv[k,4] += inc

                end

                if (k1 != 0)

                    # k1 <= 0 ; k2 >= 0
                    k = 1 + k2 + kmax + (2*kmax+1) * (-k1 + kmax)
                    # betak = -k1*beta_1 + k2*beta_2 + m*beta_3 
                    cosbetak = cosbetak1*cosbetak2*cosmbeta3 + sinbetak1*sinbetak2*cosmbeta3 + sinbetak1*cosbetak2*sinmbeta3 - cosbetak1*sinbetak2*sinmbeta3
                    inc = pref*cosbetak
                    tabWuWv[k,3] += inc * Jacv 
                    tabWuWv[k,4] += inc

                    if (k2 != 0)

                        # k1 <= 0 ; k2 <= 0
                        k = 1 - k2 + kmax + (2*kmax+1) * (-k1 + kmax)
                        # betak = -k1*beta_1 +-k2*beta_2 + m*beta_3 
                        cosbetak = cosbetak1*cosbetak2*cosmbeta3 - sinbetak1*sinbetak2*cosmbeta3 + sinbetak1*cosbetak2*sinmbeta3 + cosbetak1*sinbetak2*sinmbeta3
                        inc = pref*cosbetak
                        tabWuWv[k,3] += inc * Jacv 
                        tabWuWv[k,4] += inc
                    end
                end

            end

            cosbetak2, sinbetak2 = cosbetak2*cosbeta2 - sinbetak2*sinbeta2, sinbetak2*cosbeta2 + cosbetak2*sinbeta2

        end
        cosbetak1, sinbetak1 = cosbetak1*cosbeta1 - sinbetak1*sinbeta1, sinbetak1*cosbeta1 + cosbetak1*sinbeta1
    end

end



function tabWkp_alt_half(p::Int64, l::Int64, tabWkp::Matrix{Float64}, tabWuWv::Matrix{Float64}, m::Int64, elem::OrbitalElements_alt_half,  freq_matrix::Matrix{Float64}, nbk::Int64, nbt::Int64=nbt_default)


    duEff = pi/nbt
    dvEff = pi/2/nbt

 
    for k=1:nbk
        for i=1:4
            tabWuWv[k,i] = 0.0  

        end

    end

    alpha_1 = pi*1.0
    alpha_2 = 0.0
    alpha_3 = 0.0

    indexu = 1

     for iu=1:nbt 

        # Step 1
        Flmn_u = elem.tab_Flmn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz

        pref = 1.0/6.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_half!(m,l,tabWuWv,pref,Jacu,alpha_1,alpha_2,alpha_3)

        dalpha1_1 = duEff * dpudJu  
        dalpha1_2 = duEff * dpudJv  
        dalpha1_3 = duEff * dpudLz 



        # Step 2
        indexu += 1
        Flmn_u = elem.tab_Flmn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz

        pref = 1.0/3.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_half!(m,l,tabWuWv,pref,Jacu,alpha_1-0.5*dalpha1_1,alpha_2-0.5*dalpha1_2,alpha_3-0.5*dalpha1_3)
        
        dalpha2_1 = duEff * dpudJu  
        dalpha2_2 = duEff * dpudJv  
        dalpha2_3 = duEff * dpudLz 


        # Step 3

        pref = 1.0/3.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_half!(m,l,tabWuWv,pref,Jacu,alpha_1-0.5*dalpha2_1,alpha_2-0.5*dalpha2_2,alpha_3-0.5*dalpha2_3)
    
        dalpha3_1 = dalpha2_1
        dalpha3_2 = dalpha2_2
        dalpha3_3 = dalpha2_3

        # Step 4
        indexu += 1
        Flmn_u = elem.tab_Flmn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz


        pref = 1.0/6.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_half!(m,l,tabWuWv,pref,Jacu,alpha_1-dalpha3_1,alpha_2-dalpha3_2,alpha_3-dalpha3_3)
  
        dalpha4_1 = duEff * dpudJu  
        dalpha4_2 = duEff * dpudJv  
        dalpha4_3 = duEff * dpudLz 

        alpha_1 -= 1.0/6.0*(dalpha1_1+2.0*dalpha2_1+2.0*dalpha3_1+dalpha4_1)
        alpha_2 -= 1.0/6.0*(dalpha1_2+2.0*dalpha2_2+2.0*dalpha3_2+dalpha4_2)
        alpha_3 -= 1.0/6.0*(dalpha1_3+2.0*dalpha2_3+2.0*dalpha3_3+dalpha4_3)

        
    end

    beta_1 = 0.0
    beta_2 = pi*0.5
    beta_3 = 0.0


    indexv = 1
     for iv=1:nbt 

        # Step 1
        ylm_v = elem.tab_ylm[indexv]
        tpv = elem.tab_tpv[indexv]

        Jacv = elem.tab_Jacv[indexv]

        dpvdE = elem.tab_dpvdE[indexv]
        dpvdI3 = elem.tab_dpvdI3[indexv]
        dpvdLz = elem.tab_dpvdLz[indexv]

        dpvdJu = freq_matrix[1,1]*dpvdE + freq_matrix[2,1]*dpvdI3
        dpvdJv = freq_matrix[1,2]*dpvdE + freq_matrix[2,2]*dpvdI3
        dpvdLz = freq_matrix[1,3]*dpvdE + freq_matrix[2,3]*dpvdI3 + dpvdLz

        pref = 1.0/6.0 * 1.0/nbt * ylm_v/tpv
        increment_Wv_half!(m,l,tabWuWv,pref,Jacv,beta_1,beta_2,beta_3)


        dbeta1_1 = dvEff * dpvdJu  
        dbeta1_2 = dvEff * dpvdJv  
        dbeta1_3 = dvEff * dpvdLz 
 

        # Step 2
        indexv += 1
        ylm_v = elem.tab_ylm[indexv]
        tpv = elem.tab_tpv[indexv]

        Jacv = elem.tab_Jacv[indexv]

        dpvdE = elem.tab_dpvdE[indexv]
        dpvdI3 = elem.tab_dpvdI3[indexv]
        dpvdLz = elem.tab_dpvdLz[indexv]

        dpvdJu = freq_matrix[1,1]*dpvdE + freq_matrix[2,1]*dpvdI3
        dpvdJv = freq_matrix[1,2]*dpvdE + freq_matrix[2,2]*dpvdI3
        dpvdLz = freq_matrix[1,3]*dpvdE + freq_matrix[2,3]*dpvdI3 + dpvdLz

        pref = 1.0/3.0 * 1.0/nbt * ylm_v/tpv
        increment_Wv_half!(m,l,tabWuWv,pref,Jacv,beta_1-0.5*dbeta1_1,beta_2-0.5*dbeta1_2,beta_3-0.5*dbeta1_3)

        
        dbeta2_1 = dvEff * dpvdJu  
        dbeta2_2 = dvEff * dpvdJv  
        dbeta2_3 = dvEff * dpvdLz 

        # Step 3

        pref = 1.0/3.0 * 1.0/nbt * ylm_v/tpv
        increment_Wv_half!(m,l,tabWuWv,pref,Jacv,beta_1-0.5*dbeta2_1,beta_2-0.5*dbeta2_2,beta_3-0.5*dbeta2_3)

        dbeta3_1 = dbeta2_1 
        dbeta3_2 = dbeta2_2
        dbeta3_3 = dbeta2_3

        # Step 4
        indexv += 1
        ylm_v = elem.tab_ylm[indexv]
        tpv = elem.tab_tpv[indexv]

        Jacv = elem.tab_Jacv[indexv]

        dpvdE = elem.tab_dpvdE[indexv]
        dpvdI3 = elem.tab_dpvdI3[indexv]
        dpvdLz = elem.tab_dpvdLz[indexv]

        dpvdJu = freq_matrix[1,1]*dpvdE + freq_matrix[2,1]*dpvdI3
        dpvdJv = freq_matrix[1,2]*dpvdE + freq_matrix[2,2]*dpvdI3
        dpvdLz = freq_matrix[1,3]*dpvdE + freq_matrix[2,3]*dpvdI3 + dpvdLz

        pref = 1.0/6.0 * 1.0/nbt * ylm_v/tpv
        increment_Wv_half!(m,l,tabWuWv,pref,Jacv,beta_1-dbeta3_1,beta_2-dbeta3_2,beta_3-dbeta3_3)

        
        dbeta4_1 = dvEff * dpvdJu  
        dbeta4_2 = dvEff * dpvdJv  
        dbeta4_3 = dvEff * dpvdLz 


        beta_1 -= 1.0/6.0*(dbeta1_1+2.0*dbeta2_1+2.0*dbeta3_1+dbeta4_1)
        beta_2 -= 1.0/6.0*(dbeta1_2+2.0*dbeta2_2+2.0*dbeta3_2+dbeta4_2)
        beta_3 -= 1.0/6.0*(dbeta1_3+2.0*dbeta2_3+2.0*dbeta3_3+dbeta4_3)
       
       
    end

    # println("---")

    for k=1:nbk 

        tabWkp[p,k] = tabWuWv[k,1]*tabWuWv[k,3] + tabWuWv[k,2]*tabWuWv[k,4]


    end
    



end




# THIS IS THE CORRECT FUNCTION
# USE THIS
function ResponseMatrix_m_sampling_alt_half(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

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

    dtJu = pi/2.0/nbJ
    dtJv = pi/2.0/nbJ
    dtLz = pi/nbJ

    tab_omega = zeros(Float64, nbomega, 2) # (re, im)

    for iomega=1:nbomega 

        # iomega - 1 = im-1 + nbIm*(ire-1)

        ire = floor(Int64,(iomega-1)/nbIm) + 1
        iim = iomega - nbIm*(ire-1)

        re_omega = re_omega_min + (re_omega_max-re_omega_min)/(nbRe-1)*(ire-1)
        im_omega = im_omega_min + (im_omega_max-im_omega_min)/(nbIm-1)*(iim-1)

        tab_omega[iomega,1], tab_omega[iomega,2] = re_omega, im_omega 
    end


    elem = [OrbitalElements_alt_half_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

  

    tab_Mpq_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()] 

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 


    # Maybe use different nbJu, nbJv, nbLz
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

        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0 + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)



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

       

        for p=1:nbp 

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l 

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))
            
            

            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E,I3,Lz,u0,u1,v0,v1,grad_matrix[id_thread],nbt)
             tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)
           
           
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

                # println((omega,invden))

                for index=1:nbp*nbp 

                    # index - 1 = (q-1) + (p-1)*nbp
    
                    p = floor(Int64, (index-1)/nbp) + 1
                    q  = index  - (p-1)*nbp

    
                   

                    lp = abs(m) + floor(Int64,(p-1)/nbn)
                    lq = abs(m) + floor(Int64,(q-1)/nbn)

                    # Should use parity of l+m+k2 here ?
                    # if (mod(lp+m+k2,2) == 0 && mod(lq+m+k2,2) == 0 )

                    if (mod(lp+m+k2,2)==0 && mod(lq+m+k2,2)==0)

                        Wkp = tabWkp_temp[id_thread][p,k]  
                        Wkq = tabWkp_temp[id_thread][q,k]  
                        integrand = Wkp * Wkq * invden

                        tab_Mpq_par[index,iomega,id_thread] += integrand * jacu * jacv * jacz
                    end
 
                end

            end

        end

        # println("-----")


    end

    tab_Mpq = zeros(ComplexF64, nbp, nbp, nbomega)

 
    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz
 

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

    # pt=plot([tab_integrand[1:nbJ_default,1] tab_integrand[1:nbJ_default,2]])
    # savefig(pt, "integrand.png")

    return tab_Mpq, tab_omega
end



# THIS IS THE CORRECT FUNCTION
# USE THIS
function ResponseMatrix_m_sampling_alt_half_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default, iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default)

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

    dtJu = pi/2.0/nbJ
    dtJv = pi/2.0/nbJ
    dtLz = pi/nbJ

    tab_omega = zeros(Float64, nbomega, 2) # (re, im)

    for iomega=1:nbomega 

        # iomega - 1 = im-1 + nbIm*(ire-1)

        ire = floor(Int64,(iomega-1)/nbIm) + 1
        iim = iomega - nbIm*(ire-1)

        re_omega = re_omega_min + (re_omega_max-re_omega_min)/(nbRe-1)*(ire-1)
        im_omega = im_omega_min + (im_omega_max-im_omega_min)/(nbIm-1)*(iim-1)

        tab_omega[iomega,1], tab_omega[iomega,2] = re_omega, im_omega 
    end


    elem = [OrbitalElements_alt_half_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

  

    tab_Mpq_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()] 

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_ang_k1 = [zeros(ComplexF64, kmax+1) for it=1:Threads.nthreads()] 
    tab_ang_k2 = [zeros(ComplexF64, kmax+1) for it=1:Threads.nthreads()] 
    
    # Fill the tables of Wkp


 
    # Split the iJ threading between different nodes
    # E.g., [1,384] = [1,128] + [129,256] + [257,384] for 3 nodes
    Threads.@threads for iJ=iJmin:iJmax

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

        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0 + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)


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

        # OrbitalElements_update!(elem[id_thread],E,I3,Lz,u0,u1,v0,v1,nbt,freq_matrix[id_thread]) # l does not change: We don't have to update the spheroidal harmonics


        for p=1:nbp 

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l 

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))
            
            

            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E,I3,Lz,u0,u1,v0,v1,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)
           
           
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
    
                    integrand = Wkp * Wkq * invden

                    lp = abs(m) + floor(Int64,(p-1)/nbn)
                    lq = abs(m) + floor(Int64,(q-1)/nbn)

                    # Should use parity of l+m+k2 here ?
                    # if (mod(lp+m+k2,2) == 0 && mod(lq+m+k2,2) == 0 )

                    if (mod(lp+m+k2,2)==0 && mod(lq+m+k2,2)==0)

                            tab_Mpq_par[index,iomega,id_thread] += integrand * jacu * jacv * jacz
                    end


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
    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz
    # pref = (2*pi)^3 * 4*pi*G *dtJu*dtJv*dtLz

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


