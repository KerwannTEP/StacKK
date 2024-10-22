function increment_Wu_bare!(m::Int64, tabWuWv::Matrix{Float64}, pref::Float64, Jacu::Float64, alpha_1::Float64, alpha_2::Float64, alpha_3::Float64) #, tab_ang_k1::Array{ComplexF64}, tab_ang_k2::Array{ComplexF64})

    @fastmath sinalpha1, cosalpha1 = sincos(alpha_1)
    @fastmath sinalpha2, cosalpha2 = sincos(alpha_2)
    @fastmath sinmalpha3, cosmalpha3 = sincos(m*alpha_3)


    cosalphak1 = 1.0
    sinalphak1 = 0.0

    for k1=0:kmax 
        cosalphak2 = 1.0
        sinalphak2 = 0.0

        for k2=0:kmax

             
            
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
             

            cosalphak2, sinalphak2 = cosalphak2*cosalpha2 - sinalphak2*sinalpha2, sinalphak2*cosalpha2 + cosalphak2*sinalpha2
        end
        cosalphak1, sinalphak1 = cosalphak1*cosalpha1 - sinalphak1*sinalpha1, sinalphak1*cosalpha1 + cosalphak1*sinalpha1
 
    end

end


function increment_Wv_bare!(m::Int64, l::Int64, tabWuWv::Matrix{Float64}, pref::Float64, Jacv::Float64, beta_1::Float64, beta_2::Float64, beta_3::Float64)

    @fastmath sinbeta1, cosbeta1 = sincos(beta_1)
    @fastmath sinbeta2, cosbeta2 = sincos(beta_2)
    @fastmath sinmbeta3, cosmbeta3 = sincos(m*beta_3)

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



function tabWu_kn_bare!(n::Int64, n0::Int64, tabWu_kn::Matrix{Float64}, tabWuJ_kn::Matrix{Float64}, tabWuWv::Matrix{Float64}, m::Int64, elem::OrbitalElements_bare,  freq_matrix::Matrix{Float64}, nbk::Int64, nbt::Int64=nbt_default)


    duEff = pi/nbt
    

 
    @inbounds for k=1:nbk
        tabWuWv[k,1] = 0.0  
        tabWuWv[k,2] = 0.0  
        
    end

    alpha_1 = pi*1.0
    alpha_2 = 0.0
    alpha_3 = 0.0

    indexu = 1

    @inbounds for iu=1:nbt 

        # Step 1
        Flmn_u = elem.tab_Fn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz

        pref = 1.0/6.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_bare!(m,tabWuWv,pref,Jacu,alpha_1,alpha_2,alpha_3)

        dalpha1_1 = duEff * dpudJu  
        dalpha1_2 = duEff * dpudJv  
        dalpha1_3 = duEff * dpudLz 



        # Step 2
        indexu += 1
        Flmn_u = elem.tab_Fn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz

        pref = 1.0/3.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_bare!(m,tabWuWv,pref,Jacu,alpha_1-0.5*dalpha1_1,alpha_2-0.5*dalpha1_2,alpha_3-0.5*dalpha1_3)
        
        dalpha2_1 = duEff * dpudJu  
        dalpha2_2 = duEff * dpudJv  
        dalpha2_3 = duEff * dpudLz 


        # Step 3

        pref = 1.0/3.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_bare!(m,tabWuWv,pref,Jacu,alpha_1-0.5*dalpha2_1,alpha_2-0.5*dalpha2_2,alpha_3-0.5*dalpha2_3)
    
        dalpha3_1 = dalpha2_1
        dalpha3_2 = dalpha2_2
        dalpha3_3 = dalpha2_3

        # Step 4
        indexu += 1
        Flmn_u = elem.tab_Fn[indexu]
        tpu = elem.tab_tpu[indexu]
        Jacu = elem.tab_Jacu[indexu]

        dpudE = elem.tab_dpudE[indexu]
        dpudI3 = elem.tab_dpudI3[indexu]
        dpudLz = elem.tab_dpudLz[indexu]

        dpudJu = freq_matrix[1,1]*dpudE + freq_matrix[2,1]*dpudI3
        dpudJv = freq_matrix[1,2]*dpudE + freq_matrix[2,2]*dpudI3
        dpudLz = freq_matrix[1,3]*dpudE + freq_matrix[2,3]*dpudI3 + dpudLz


        pref = 1.0/6.0 * 1.0/nbt * Flmn_u/tpu
        increment_Wu_bare!(m,tabWuWv,pref,Jacu,alpha_1-dalpha3_1,alpha_2-dalpha3_2,alpha_3-dalpha3_3)
  
        dalpha4_1 = duEff * dpudJu  
        dalpha4_2 = duEff * dpudJv  
        dalpha4_3 = duEff * dpudLz 

        alpha_1 -= 1.0/6.0*(dalpha1_1+2.0*dalpha2_1+2.0*dalpha3_1+dalpha4_1)
        alpha_2 -= 1.0/6.0*(dalpha1_2+2.0*dalpha2_2+2.0*dalpha3_2+dalpha4_2)
        alpha_3 -= 1.0/6.0*(dalpha1_3+2.0*dalpha2_3+2.0*dalpha3_3+dalpha4_3)

        
    end

    @inbounds for k=1:nbk 

        tabWu_kn[n-n0+1,k] = tabWuWv[k,1]
        tabWuJ_kn[n-n0+1,k] = tabWuWv[k,2]
        


    end

    

end

function tabWv_kl_bare!(l::Int64, tabWv_kl::Matrix{Float64}, tabWvJ_kl::Matrix{Float64}, tabWuWv::Matrix{Float64}, m::Int64, elem::OrbitalElements_bare,  freq_matrix::Matrix{Float64}, nbk::Int64, nbt::Int64=nbt_default)


    dvEff = pi/2/nbt

    beta_1 = 0.0
    beta_2 = pi*0.5
    beta_3 = 0.0

    @inbounds for k=1:nbk
        tabWuWv[k,3] = 0.0  
        tabWuWv[k,4] = 0.0  

    end


    indexv = 1
    @inbounds for iv=1:nbt 

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
        increment_Wv_bare!(m,l,tabWuWv,pref,Jacv,beta_1,beta_2,beta_3)


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
        increment_Wv_bare!(m,l,tabWuWv,pref,Jacv,beta_1-0.5*dbeta1_1,beta_2-0.5*dbeta1_2,beta_3-0.5*dbeta1_3)

        
        dbeta2_1 = dvEff * dpvdJu  
        dbeta2_2 = dvEff * dpvdJv  
        dbeta2_3 = dvEff * dpvdLz 

        # Step 3

        pref = 1.0/3.0 * 1.0/nbt * ylm_v/tpv
        increment_Wv_bare!(m,l,tabWuWv,pref,Jacv,beta_1-0.5*dbeta2_1,beta_2-0.5*dbeta2_2,beta_3-0.5*dbeta2_3)

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
        increment_Wv_bare!(m,l,tabWuWv,pref,Jacv,beta_1-dbeta3_1,beta_2-dbeta3_2,beta_3-dbeta3_3)

        
        dbeta4_1 = dvEff * dpvdJu  
        dbeta4_2 = dvEff * dpvdJv  
        dbeta4_3 = dvEff * dpvdLz 


        beta_1 -= 1.0/6.0*(dbeta1_1+2.0*dbeta2_1+2.0*dbeta3_1+dbeta4_1)
        beta_2 -= 1.0/6.0*(dbeta1_2+2.0*dbeta2_2+2.0*dbeta3_2+dbeta4_2)
        beta_3 -= 1.0/6.0*(dbeta1_3+2.0*dbeta2_3+2.0*dbeta3_3+dbeta4_3)
       
       
    end

    # println("---")

    @inbounds for k=1:nbk 

        tabWvJ_kl[l-abs(m)+1,k] = tabWuWv[k,3]
        tabWv_kl[l-abs(m)+1,k] = tabWuWv[k,4]
        


    end
    



end


function tab_Mpq_bare_par!(id_thread::Int64, tab_Mpq_par::Matrix{ComplexF64}, tabWu_kn_temp::Matrix{Float64}, tabWuJ_kn_temp::Matrix{Float64}, tabWv_kl_temp::Matrix{Float64}, tabWvJ_kl_temp::Matrix{Float64}, omega::ComplexF64, k::Int64, m::Int64, dFdJu::Float64, dFdJv::Float64, dFdLz::Float64, Omegau::Float64, Omegav::Float64, Omegaz::Float64, jacu::Float64, jacv::Float64, jacz::Float64, nbn::Int64, nbp::Int64, n0::Int64)

    # nbkline = 2*kmax+1
    # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)
    
    k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
    k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

    kdotdF = k1*dFdJu + k2*dFdJv + m*dFdLz
    kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz 

    invden = kdotdF/(omega-kdotOmega) * jacu * jacv * jacz


    for p=1:nbp 

        lp = abs(m) + floor(Int64,(p-1)/nbn)
        np = p - 1 - (lp-abs(m))*nbn + n0
        Wkp = tabWu_kn_temp[np-n0+1,k]*tabWvJ_kl_temp[lp-abs(m)+1,k] + tabWuJ_kn_temp[np-n0+1,k]*tabWv_kl_temp[lp-abs(m)+1,k] 

        for q=1:p-1

            lq = abs(m) + floor(Int64,(q-1)/nbn)
            
            if (mod(lp+m+k2,2)==0 && mod(lq+m+k2,2)==0)

                nq = q - 1 - (lq-abs(m))*nbn + n0
                Wkq = tabWu_kn_temp[nq-n0+1,k]*tabWvJ_kl_temp[lq-abs(m)+1,k] + tabWuJ_kn_temp[nq-n0+1,k]*tabWv_kl_temp[lq-abs(m)+1,k] 
                index = 1 + (q-1) + (p-1)*nbp
                indexSym = 1 + (p-1) + (q-1)*nbp

                integrand = Wkp * Wkq * invden

                tab_Mpq_par[index,id_thread] += integrand 
                tab_Mpq_par[indexSym,id_thread] += integrand 


            end

        end

        # case p=q

        if (mod(lp+m+k2,2)==0)

            index = 1 + (p-1) + (p-1)*nbp
            integrand = Wkp * Wkp * invden 

            tab_Mpq_par[index,id_thread] += integrand 

        end
    end

    

end


function fill_W!(elem::OrbitalElements_bare, E::Float64, I3::Float64, Lz::Float64, u0::Float64, u1::Float64,
                grad_matrix::Matrix{Float64}, nbt::Int64, n0::Int64, tabWu_kn_temp::Matrix{Float64},
                tabWuJ_kn_temp::Matrix{Float64}, tabWuWv::Matrix{Float64}, m::Int64, freq_matrix::Matrix{Float64},
                nbk::Int64, v0::Float64, v1::Float64, tabWv_kl_temp::Matrix{Float64}, tabWvJ_kl_temp::Matrix{Float64})

    for n=n0:nmax 

        OrbitalElements_bare_update_u!(elem,n,E,I3,Lz,u0,u1,grad_matrix,nbt)
        tabWu_kn_bare!(n,n0,tabWu_kn_temp,tabWuJ_kn_temp,tabWuWv,m,elem,freq_matrix,nbk,nbt)

    end

    for l=abs(m):lmax 

        OrbitalElements_bare_update_v!(elem,l,m,E,I3,Lz,v0,v1,grad_matrix,nbt)
        tabWv_kl_bare!(l,tabWv_kl_temp,tabWvJ_kl_temp,tabWuWv,m,elem,freq_matrix,nbk,nbt)

    end

end

function tab_Mpq_bare!(index::Int64, tab_Mpq_bare::Matrix{ComplexF64}, tab_Mpq_bare_par::Matrix{ComplexF64}, nbp::Int64, pref::Float64)

    # index - 1 = (q-1) + (p-1)*nbp

    p = floor(Int64, (index-1)/nbp) + 1
    q  = index  - (p-1)*nbp


    for it=1:Threads.nthreads()

        tab_Mpq_bare[p,q] += tab_Mpq_bare_par[index,it]

    end
    tab_Mpq_bare[p,q] *= pref

end



function ResponseMatrix_m_bare(m::Int64, omega::ComplexF64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

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

    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJ
    dtJv = pi/2.0*(1.0-epsJ)/nbJ
    dtLz = pi*(1.0-epsJ)/nbJ

    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

    tab_Mpq_bare_par = [0.0 + 0.0im for index=1:nbp*nbp, it=1:Threads.nthreads()] 

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()] 

    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

 
    # Fill the tables of Wkp

    Threads.@threads for iJ=1:nbgrid

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ 
        # nbv = nbJ 

        # println("----")

        iu = floor(Int64,(iJ - 1)/nbJ^2) + 1
        leftover = iJ - 1 - nbJ^2*(iu-1)

        # leftover = iz-1 + nbz*(iv-1)
        # nbz = nbJ 

        iv = floor(Int64,leftover/nbJ) + 1
        iz = leftover - nbJ*(iv-1) + 1

        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0+epsJ + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)

        # println("----")

        E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

        u0, u1 = find_bounds_u(E,Lz,I3)
        v0, v1 = find_bounds_v(E,Lz,I3)

        id_thread = Threads.threadid()

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E,Lz,I3)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  


        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        fill_W!(elem[id_thread],E,I3,Lz,u0,u1,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0,v1,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])


        # invert k and pq sum ?
        for k=1:nbk

            tab_Mpq_bare_par!(id_thread,tab_Mpq_bare_par,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],omega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0)

        end

    end

    

    tab_Mpq_bare = zeros(ComplexF64, nbp, nbp)

    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare!(index,tab_Mpq_bare,tab_Mpq_bare_par,nbp,pref)
       
    end

    return tab_Mpq_bare
end




function ResponseMatrix_m_bare_separate(m::Int64, omega::ComplexF64, nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=25, nbt::Int64=nbt_default)

    @assert (abs(m) <= mmax) "m must be less that m_max"

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbgrid = nbJu*nbJv*nbLz
    nbk = (2*kmax+1)^2

    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJu
    dtJv = pi/2.0*(1.0-epsJ)/nbJv
    dtLz = pi*(1.0-epsJ)/nbLz

    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

    tab_Mpq_bare_par = [0.0 + 0.0im for index=1:nbp*nbp, it=1:Threads.nthreads()] 

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()] 

    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

 
    # Fill the tables of Wkp

    Threads.@threads for iJ=1:nbgrid

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ 
        # nbv = nbJ 

        # println("----")

        iu = floor(Int64,(iJ - 1)/(nbJv*nbLz)) + 1
        leftover = iJ - 1 - nbJv*nbLz*(iu-1)

        # leftover = iz-1 + nbz*(iv-1)
        # nbz = nbJ 

        iv = floor(Int64,leftover/nbLz) + 1
        iz = leftover - nbLz*(iv-1) + 1

        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0+epsJ + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)

        # println("----")

        E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

        u0, u1 = find_bounds_u(E,Lz,I3)
        v0, v1 = find_bounds_v(E,Lz,I3)

        id_thread = Threads.threadid()

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E,Lz,I3)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  


        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        fill_W!(elem[id_thread],E,I3,Lz,u0,u1,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0,v1,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])


        # invert k and pq sum ?
        for k=1:nbk

            tab_Mpq_bare_par!(id_thread,tab_Mpq_bare_par,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],omega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0)

        end

    end

    

    tab_Mpq_bare = zeros(ComplexF64, nbp, nbp)

    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare!(index,tab_Mpq_bare,tab_Mpq_bare_par,nbp,pref)
       
    end

    return tab_Mpq_bare
end



function compute_tab_Mpq!(m::Int64, n0::Int64, index::Int64, nbn::Int64, nbp::Int64, 
                        tab_Mpq::Matrix{ComplexF64}, tab_Mpq_bare::Matrix{ComplexF64})

    
    # index - 1 = (q-1) + (p-1)*nbp                    
    
    p = floor(Int64, (index-1)/nbp) + 1
    q  = index  - (p-1)*nbp

    lp = abs(m) + floor(Int64,(p-1)/nbn)
    lq = abs(m) + floor(Int64,(q-1)/nbn)

    np = p - 1 - (lp-abs(m))*nbn + n0
    nq = q - 1 - (lq-abs(m))*nbn + n0

    id_p = tab_index[lp+1,abs(m)+1]
    id_q = tab_index[lq+1,abs(m)+1]

    for i=n0:nmax

        p_i = 1 + (lp - abs(m))*nbn + i-n0
        coefGS_i = tab_YGS[id_p][i-n0+1,np-n0+1]

        for j=n0:nmax

            q_j = 1 + (lq - abs(m))*nbn + j-n0
            coefGS_j = tab_YGS[id_q][j-n0+1,nq-n0+1]

            M_pi_qj_bare = tab_Mpq_bare[p_i,q_j]
        
            # Delta because of normalization of bi-orthogonal basis elements
            tab_Mpq[p,q] += coefGS_i*coefGS_j*M_pi_qj_bare*Delta
        end

    end

end





function ResponseMatrix_m_GS(m::Int64, omega::ComplexF64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn

    tab_Mpq_bare = ResponseMatrix_m_bare(m,omega,nbJ,nbt)


    tab_Mpq = zeros(ComplexF64, nbp, nbp)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq!(m,n0,index,nbn,nbp,tab_Mpq,tab_Mpq_bare)

    end

    return tab_Mpq

end

function ResponseMatrix_m_GS_separate(m::Int64, omega::ComplexF64, nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, nbt::Int64=nbt_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn

    tab_Mpq_bare = ResponseMatrix_m_bare_separate(m,omega,nbJu,nbJv,nbLz,nbt)


    tab_Mpq = zeros(ComplexF64, nbp, nbp)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq!(m,n0,index,nbn,nbp,tab_Mpq,tab_Mpq_bare)

    end

    return tab_Mpq

end

############################################################################################
# Sample the response matrix over a grid of complex frequencies
############################################################################################

function tab_Mpq_bare_par_sampling!(tab_Mpq_par::Matrix{ComplexF64}, tab_W_temp::Vector{Float64}, tabWu_kn_temp::Matrix{Float64}, tabWuJ_kn_temp::Matrix{Float64}, tabWv_kl_temp::Matrix{Float64}, tabWvJ_kl_temp::Matrix{Float64}, 
                                    tab_omega::Matrix{Float64}, nbomega::Int64, k::Int64, m::Int64, dFdJu::Float64, dFdJv::Float64, dFdLz::Float64, Omegau::Float64, Omegav::Float64, Omegaz::Float64, jacu::Float64, 
                                    jacv::Float64, jacz::Float64, nbn::Int64, nbp::Int64, n0::Int64)

    # nbkline = 2*kmax+1
    # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

    
    
    k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
    k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

    kdotdF = k1*dFdJu + k2*dFdJv + m*dFdLz
    kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz 


    for p=1:nbp 

        lp = abs(m) + floor(Int64,(p-1)/nbn)
        if (mod(lp+m+k2,2)==0)
            np = p - 1 - (lp-abs(m))*nbn + n0
            @inbounds Wkp = tabWu_kn_temp[np-n0+1,k]*tabWvJ_kl_temp[lp-abs(m)+1,k] + tabWuJ_kn_temp[np-n0+1,k]*tabWv_kl_temp[lp-abs(m)+1,k] 
            tab_W_temp[p] = Wkp 
        end
    end

 
    for p=1:nbp 

        lp = abs(m) + floor(Int64,(p-1)/nbn)

        if (mod(lp+m+k2,2)==0)

            # np = p - 1 - (lp-abs(m))*nbn + n0
            
            Wkp = tab_W_temp[p]

            pref = Wkp * kdotdF * jacu * jacv * jacz 


            for q=1:p-1

                lq = abs(m) + floor(Int64,(q-1)/nbn)
                
                if (mod(lq+m+k2,2)==0)

                    # nq = q - 1 - (lq-abs(m))*nbn + n0
                    Wkq = tab_W_temp[q]
                    
                    
                    index = 1 + (q-1) + (p-1)*nbp
                    indexSym = 1 + (p-1) + (q-1)*nbp


                    for iomega=1:nbomega

                        @inbounds re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                        omega = re_omega + 1.0im*im_omega
                
                        den = (omega-kdotOmega) 
                        integrand = pref * Wkq / den

                        @inbounds tab_Mpq_par[iomega,index] += integrand 
                        @inbounds tab_Mpq_par[iomega,indexSym] += integrand 
                    end

                end

            end

            # case p=q

            index = 1 + (p-1) + (p-1)*nbp     

            for iomega=1:nbomega

                @inbounds re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
        
                den = (omega-kdotOmega) 
                integrand = Wkp * Wkp * kdotdF * jacu * jacv * jacz / den

                @inbounds tab_Mpq_par[iomega,index] += integrand 
            end


        end

    end


end

function tab_Mpq_bare_sampling!(index::Int64, tab_Mpq_bare::Array{ComplexF64, 3}, tab_Mpq_bare_par::Vector{Matrix{ComplexF64}}, nbp::Int64, nbomega::Int64, pref::Float64)

    # index - 1 = (q-1) + (p-1)*nbp

    p = floor(Int64, (index-1)/nbp) + 1
    q  = index  - (p-1)*nbp

    for iomega=1:nbomega

        for it=1:Threads.nthreads()
            @inbounds tab_Mpq_bare[p,q,iomega] += tab_Mpq_bare_par[it][iomega,index]

        end
        @inbounds tab_Mpq_bare[p,q,iomega] *= pref

    end


end


function ResponseMatrix_m_bare_sampling(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
                                        nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

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

    tab_omega = zeros(Float64, nbomega, 2) # (re, im)

    for iomega=1:nbomega 

        # iomega - 1 = im-1 + nbIm*(ire-1)

        ire = floor(Int64,(iomega-1)/nbIm) + 1
        iim = iomega - nbIm*(ire-1)

        if (nbRe > 1)
            re_omega = re_omega_min + (re_omega_max-re_omega_min)/(nbRe-1)*(ire-1)
        else
            re_omega = re_omega_min
        end

        if (nbIm > 1)
            im_omega = im_omega_min + (im_omega_max-im_omega_min)/(nbIm-1)*(iim-1)
        else
            im_omega = im_omega_min
        end

        tab_omega[iomega,1], tab_omega[iomega,2] = re_omega, im_omega 
    end


    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJ
    dtJv = pi/2.0*(1.0-epsJ)/nbJ
    dtLz = pi*(1.0-epsJ)/nbJ

    # dtJu = Jumax/nbJ
    # dtJv = Jvmax/nbJ
    # dtLz = 2.0*Lzmax/nbJ





    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

    # tab_Mpq_bare_par = [zeros(ComplexF64, nbp*nbp, nbomega) for it=1:Threads.nthreads()] 
    tab_Mpq_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 

    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_W_temp = [zeros(Float64, nbp) for it=1:Threads.nthreads()] 
 
    # Fill the tables of Wkp

    Threads.@threads for iJ=1:nbgrid

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ 
        # nbv = nbJ 

        # println("----")

        iu = floor(Int64,(iJ - 1)/nbJ^2) + 1
        leftover = iJ - 1 - nbJ^2*(iu-1)

        # leftover = iz-1 + nbz*(iv-1)
        # nbz = nbJ 

        iv = floor(Int64,leftover/nbJ) + 1
        iz = leftover - nbJ*(iv-1) + 1





        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0 + epsJ + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)


        # Ju = dtJu*(iu-0.5)
        # Jv = dtJv*(iv-0.5)
        # Lz = -Lzmax + dtLz*(iz-0.5)

        # jacu = 1.0
        # jacv = 1.0
        # jacz = 1.0





        # println("----")

        E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

        u0, u1 = find_bounds_u(E,Lz,I3)
        v0, v1 = find_bounds_v(E,Lz,I3)

        id_thread = Threads.threadid()

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E,Lz,I3)
        @inbounds Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  


        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        fill_W!(elem[id_thread],E,I3,Lz,u0,u1,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0,v1,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])


        for k=1:nbk
            tab_Mpq_bare_par_sampling!(tab_Mpq_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0)

        end

    end

    tab_Mpq_bare = zeros(ComplexF64, nbp, nbp, nbomega)

    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare_sampling!(index,tab_Mpq_bare,tab_Mpq_bare_par,nbp,nbomega,pref)
       
    end

    return tab_Mpq_bare, tab_omega
end

function compute_tab_Mpq_sampling!(m::Int64, n0::Int64, index::Int64, nbn::Int64, nbp::Int64, 
                        tab_Mpq::Array{ComplexF64}, tab_Mpq_bare::Array{ComplexF64}, nbomega::Int64)


    # index - 1 = (q-1) + (p-1)*nbp                    

    p = floor(Int64, (index-1)/nbp) + 1
    q  = index  - (p-1)*nbp

    lp = abs(m) + floor(Int64,(p-1)/nbn)
    lq = abs(m) + floor(Int64,(q-1)/nbn)

    np = p - 1 - (lp-abs(m))*nbn + n0
    nq = q - 1 - (lq-abs(m))*nbn + n0

    id_p = tab_index[lp+1,abs(m)+1]
    id_q = tab_index[lq+1,abs(m)+1]

    for i=n0:nmax

        p_i = 1 + (lp - abs(m))*nbn + i-n0
        coefGS_i = tab_YGS[id_p][i-n0+1,np-n0+1]

        for j=n0:nmax

            q_j = 1 + (lq - abs(m))*nbn + j-n0
            coefGS_j = tab_YGS[id_q][j-n0+1,nq-n0+1]

            for iomega=1:nbomega

                M_pi_qj_bare = tab_Mpq_bare[p_i,q_j,iomega]

                # Delta because of normalization of bi-orthogonal basis elements
                tab_Mpq[p,q,iomega] += coefGS_i*coefGS_j*M_pi_qj_bare*Delta

            end
        end

    end

end


function ResponseMatrix_m_GS_sampling(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
                            nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbomega = nbRe * nbIm 

    tab_Mpq_bare, tab_omega = ResponseMatrix_m_bare_sampling(m,re_omega_min,re_omega_max,im_omega_min,im_omega_max,nbRe,nbIm,nbJ,nbt)


    tab_Mpq = zeros(ComplexF64, nbp, nbp, nbomega)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq_sampling!(m,n0,index,nbn,nbp,tab_Mpq,tab_Mpq_bare,nbomega)

    end

    return tab_Mpq, tab_omega

end





############################################################################################
# Split the action integral to distribute the load over multiple nodes
############################################################################################

function ResponseMatrix_m_bare_sampling_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

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

    tab_omega = zeros(Float64, nbomega, 2) # (re, im)

    for iomega=1:nbomega 

    # iomega - 1 = im-1 + nbIm*(ire-1)

        ire = floor(Int64,(iomega-1)/nbIm) + 1
        iim = iomega - nbIm*(ire-1)

        if (nbRe > 1)
            re_omega = re_omega_min + (re_omega_max-re_omega_min)/(nbRe-1)*(ire-1)
        else
            re_omega = re_omega_min
        end

        if (nbIm > 1)
            im_omega = im_omega_min + (im_omega_max-im_omega_min)/(nbIm-1)*(iim-1)
        else
            im_omega = im_omega_min
        end

        tab_omega[iomega,1], tab_omega[iomega,2] = re_omega, im_omega 
    end


    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJ
    dtJv = pi/2.0*(1.0-epsJ)/nbJ
    dtLz = pi*(1.0-epsJ)/nbJ

    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]

    # tab_Mpq_bare_par = [zeros(ComplexF64, nbp*nbp, nbomega) for it=1:Threads.nthreads()] 
    tab_Mpq_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 

    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_W_temp = [zeros(Float64, nbp) for it=1:Threads.nthreads()] 

    # Fill the tables of Wkp

    Threads.@threads for iJ=iJmin_default:iJmax_default

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ 
        # nbv = nbJ 

        # println("----")

        iu = floor(Int64,(iJ - 1)/nbJ^2) + 1
        leftover = iJ - 1 - nbJ^2*(iu-1)

        # leftover = iz-1 + nbz*(iv-1)
        # nbz = nbJ 

        iv = floor(Int64,leftover/nbJ) + 1
        iz = leftover - nbJ*(iv-1) + 1

        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)
        tLz = -pi/2.0 + epsJ + dtLz*(iz-0.5)

        Ju = tan(tJu)
        Jv = tan(tJv)
        Lz = tan(tLz)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)
        jacz = (1.0+Lz^2)

        # println("----")

        E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

        u0, u1 = find_bounds_u(E,Lz,I3)
        v0, v1 = find_bounds_v(E,Lz,I3)

        id_thread = Threads.threadid()

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E,Lz,I3)
        @inbounds Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  


        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        fill_W!(elem[id_thread],E,I3,Lz,u0,u1,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0,v1,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])


        for k=1:nbk
            tab_Mpq_bare_par_sampling!(tab_Mpq_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0)

        end

    end

    tab_Mpq_bare = zeros(ComplexF64, nbp, nbp, nbomega)

    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare_sampling!(index,tab_Mpq_bare,tab_Mpq_bare_par,nbp,nbomega,pref)

    end

    return tab_Mpq_bare, tab_omega
end


function ResponseMatrix_m_GS_sampling_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbomega = nbRe * nbIm 

    tab_Mpq_bare, tab_omega = ResponseMatrix_m_bare_sampling_split(m,re_omega_min,re_omega_max,im_omega_min,im_omega_max,nbRe,nbIm,nbJ,nbt)


    tab_Mpq = zeros(ComplexF64, nbp, nbp, nbomega)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq_sampling!(m,n0,index,nbn,nbp,tab_Mpq,tab_Mpq_bare,nbomega)

    end

    return tab_Mpq, tab_omega

end
