function tab_Mpq_ab_bare_par_sampling_rot!(tab_Mpq_a_par::Matrix{ComplexF64}, tab_Mpq_b_par::Matrix{ComplexF64}, tab_W_temp::Vector{Float64}, tabWu_kn_temp::Matrix{Float64}, tabWuJ_kn_temp::Matrix{Float64}, tabWv_kl_temp::Matrix{Float64}, tabWvJ_kl_temp::Matrix{Float64}, 
    tab_omega::Matrix{Float64}, nbomega::Int64, k::Int64, m::Int64, dFdJu::Float64, dFdJv::Float64, dFdLz::Float64, Omegau::Float64, Omegav::Float64, Omegaz::Float64, jacu::Float64, 
    jacv::Float64, jacz::Float64, nbn::Int64, nbp::Int64, n0::Int64, Lz::Float64)

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

                        @inbounds tab_Mpq_a_par[iomega,index] += integrand 
                        @inbounds tab_Mpq_a_par[iomega,indexSym] += integrand 

                        @inbounds tab_Mpq_b_par[iomega,index] += integrand * sign(Lz)
                        @inbounds tab_Mpq_b_par[iomega,indexSym] += integrand * sign(Lz)
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

                @inbounds tab_Mpq_a_par[iomega,index] += integrand 

                @inbounds tab_Mpq_b_par[iomega,index] += integrand * sign(Lz)
            end


        end

    end


end


function tab_Mpq_c_bare_par_sampling_rot!(tab_Mpq_c_par::Matrix{ComplexF64}, tab_W_temp::Vector{Float64}, tabWu_kn_temp::Matrix{Float64}, tabWuJ_kn_temp::Matrix{Float64}, tabWv_kl_temp::Matrix{Float64}, tabWvJ_kl_temp::Matrix{Float64}, 
                                    tab_omega::Matrix{Float64}, nbomega::Int64, k::Int64, m::Int64, Ftot::Float64, Omegau::Float64, Omegav::Float64, Omegaz::Float64, jacu::Float64, 
                                    jacv::Float64, nbn::Int64, nbp::Int64, n0::Int64)

    # nbkline = 2*kmax+1
    # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

    
    
    k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
    k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

    mF = m*Ftot

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

            pref = Wkp * mF * jacu * jacv


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

                        @inbounds tab_Mpq_c_par[iomega,index] += 0.5 * integrand 
                        @inbounds tab_Mpq_c_par[iomega,indexSym] += 0.5 * integrand 
                    end

                end

            end

            # case p=q

            index = 1 + (p-1) + (p-1)*nbp     

            for iomega=1:nbomega

                @inbounds re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
        
                den = (omega-kdotOmega) 
                integrand = Wkp * Wkp * mF * jacu * jacv / den

                @inbounds tab_Mpq_c_par[iomega,index] += 0.5 * integrand 
            end


        end

    end


end

function ResponseMatrix_m_bare_sampling_rot_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

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


    dtJu = pi/2.0/nbJ
    dtJv = pi/2.0/nbJ
    dtLz = pi/nbJ

    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]


    tab_Mpq_a_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 
    tab_Mpq_b_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 


    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_W_temp = [zeros(Float64, nbp) for it=1:Threads.nthreads()] 


    Threads.@threads for iJ=iJmin:iJmax

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
        tLz = -pi/2.0 + dtLz*(iz-0.5)

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

        # implement rotation after this comment

        for k=1:nbk
            tab_Mpq_ab_bare_par_sampling_rot!(tab_Mpq_a_bare_par[id_thread],tab_Mpq_b_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0,Lz)

        end

    end


    tab_Mpq_c_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 


    Threads.@threads for iJ=iJmin2d:iJmax2d

        # iJ - 1 = (iv-1) + nbv*(iu-1)

        # TODO
        iu = floor(Int64,(iJ-1)/nbJ) + 1
        iv = iJ - 1 - nbJ*(iu-1) + 1


        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)


        Lz_p = epsLz
        Lz_m = -epsLz

        Ju = tan(tJu)
        Jv = tan(tJv)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)


        E_p, _, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz_p,Jv)
        E_m, _, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz_m,Jv)

        u0_p, u1_p = find_bounds_u(E_p,Lz_p,I3_p)
        v0_p, v1_p = find_bounds_v(E_p,Lz_p,I3_p)

        u0_m, u1_m = find_bounds_u(E_m,Lz_m,I3_m)
        v0_m, v1_m = find_bounds_v(E_m,Lz_m,I3_m)

        id_thread = Threads.threadid()

        # Lz = 0+

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_p,Lz_p,I3_p)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_p,Lz_p)

        fill_W!(elem[id_thread],E_p,I3_p,Lz_p,u0_p,u1_p,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0_p,v1_p,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])

        for k=1:nbk

            tab_Mpq_c_bare_par_sampling_rot!(tab_Mpq_c_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,Ftot,Omegau,Omegav,Omegaz,jacu,jacv,nbn,nbp,n0)

        end

        # Lz = 0-

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_m,Lz_m,I3_m)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_m,Lz_m)

        fill_W!(elem[id_thread],E_m,I3_m,Lz_m,u0_m,u1_m,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0_m,v1_m,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])

        for k=1:nbk

            tab_Mpq_c_bare_par_sampling_rot!(tab_Mpq_c_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,Ftot,Omegau,Omegav,Omegaz,jacu,jacv,nbn,nbp,n0)

        end
    end

    tab_Mpq_a_bare = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b_bare = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c_bare = zeros(ComplexF64, nbp, nbp, nbomega)

    pref  = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz
    prefc = 2.0*(2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare_sampling!(index,tab_Mpq_a_bare,tab_Mpq_a_bare_par,nbp,nbomega,pref)
        tab_Mpq_bare_sampling!(index,tab_Mpq_b_bare,tab_Mpq_b_bare_par,nbp,nbomega,pref)
        tab_Mpq_bare_sampling!(index,tab_Mpq_c_bare,tab_Mpq_c_bare_par,nbp,nbomega,prefc)

    end

    return tab_Mpq_a_bare, tab_Mpq_b_bare, tab_Mpq_c_bare, tab_omega
end

function compute_tab_Mpq_sampling_rot!(m::Int64, n0::Int64, index::Int64, nbn::Int64, nbp::Int64, 
    tab_Mpq_a::Array{ComplexF64}, tab_Mpq_b::Array{ComplexF64}, tab_Mpq_c::Array{ComplexF64}, 
    tab_Mpq_a_bare::Array{ComplexF64}, tab_Mpq_b_bare::Array{ComplexF64}, tab_Mpq_c_bare::Array{ComplexF64}, nbomega::Int64)


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

                M_a_pi_qj_bare = tab_Mpq_a_bare[p_i,q_j,iomega]
                M_b_pi_qj_bare = tab_Mpq_b_bare[p_i,q_j,iomega]
                M_c_pi_qj_bare = tab_Mpq_c_bare[p_i,q_j,iomega]

                # Delta because of normalization of bi-orthogonal basis elements
                tab_Mpq_a[p,q,iomega] += coefGS_i*coefGS_j*M_a_pi_qj_bare*Delta
                tab_Mpq_b[p,q,iomega] += coefGS_i*coefGS_j*M_b_pi_qj_bare*Delta
                tab_Mpq_c[p,q,iomega] += coefGS_i*coefGS_j*M_c_pi_qj_bare*Delta

            end
        end

    end

end

function ResponseMatrix_m_GS_sampling_rot_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbomega = nbRe * nbIm 

    tab_Mpq_a_bare, tab_Mpq_b_bare, tab_Mpq_c_bare, tab_omega = ResponseMatrix_m_bare_sampling_rot_split(m,re_omega_min,re_omega_max,im_omega_min,im_omega_max,nbRe,nbIm,nbJ,nbt,epsLz,iJmin,iJmax,iJmin2d,iJmax2d)

    tab_Mpq_a = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c = zeros(ComplexF64, nbp, nbp, nbomega)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq_sampling_rot!(m,n0,index,nbn,nbp,tab_Mpq_a,tab_Mpq_b,tab_Mpq_c,tab_Mpq_a_bare,tab_Mpq_b_bare,tab_Mpq_c_bare,nbomega)

    end

    return tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega

end


#################
#################


function ResponseMatrix_m_bare_sampling_rot_split_separate(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64, nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, 
    nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

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

    dtJu = pi/2.0*(1.0-epsJ)/nbJu
    dtJv = pi/2.0*(1.0-epsJ)/nbJv
    dtLz = pi*(1.0-epsJ)/nbLz

    elem = [OrbitalElements_bare_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]


    tab_Mpq_a_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 
    tab_Mpq_b_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 


    tabWu_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 
    tabWuJ_kn_temp = [zeros(Float64, nbn, nbk) for it=1:Threads.nthreads()] 

    tabWv_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 
    tabWvJ_kl_temp = [zeros(Float64, nbl, nbk) for it=1:Threads.nthreads()] 

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()] 

    tab_W_temp = [zeros(Float64, nbp) for it=1:Threads.nthreads()] 


    Threads.@threads for iJ=iJmin:iJmax

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
        @inbounds Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]  


        dFdE = _dFdE(E,Lz)
        dFdLz = _dFdLz(E,Lz)

        dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
        dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
        dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

        fill_W!(elem[id_thread],E,I3,Lz,u0,u1,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0,v1,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])

        # implement rotation after this comment

        for k=1:nbk
            tab_Mpq_ab_bare_par_sampling_rot!(tab_Mpq_a_bare_par[id_thread],tab_Mpq_b_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,dFdJu,dFdJv,dFdLz,Omegau,Omegav,Omegaz,jacu,jacv,jacz,nbn,nbp,n0,Lz)

        end

    end


    tab_Mpq_c_bare_par = [zeros(ComplexF64, nbomega, nbp*nbp) for it=1:Threads.nthreads()] 


    Threads.@threads for iJ=iJmin2d:iJmax2d

        # iJ - 1 = (iv-1) + nbv*(iu-1)

        # TODO
        iu = floor(Int64,(iJ-1)/nbJv) + 1
        iv = iJ - 1 - nbJv*(iu-1) + 1


        tJu = dtJu*(iu-0.5)
        tJv = dtJv*(iv-0.5)


        Lz_p = epsLz
        Lz_m = -epsLz

        Ju = tan(tJu)
        Jv = tan(tJv)

        jacu = (1.0+Ju^2)
        jacv = (1.0+Jv^2)


        E_p, _, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz_p,Jv)
        E_m, _, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz_m,Jv)

        u0_p, u1_p = find_bounds_u(E_p,Lz_p,I3_p)
        v0_p, v1_p = find_bounds_v(E_p,Lz_p,I3_p)

        u0_m, u1_m = find_bounds_u(E_m,Lz_m,I3_m)
        v0_m, v1_m = find_bounds_v(E_m,Lz_m,I3_m)

        id_thread = Threads.threadid()

        # Lz = 0+

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_p,Lz_p,I3_p)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_p,Lz_p)

        fill_W!(elem[id_thread],E_p,I3_p,Lz_p,u0_p,u1_p,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0_p,v1_p,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])

        for k=1:nbk

            tab_Mpq_c_bare_par_sampling_rot!(tab_Mpq_c_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,Ftot,Omegau,Omegav,Omegaz,jacu,jacv,nbn,nbp,n0)

        end

        # Lz = 0-

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_m,Lz_m,I3_m)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_m,Lz_m)

        fill_W!(elem[id_thread],E_m,I3_m,Lz_m,u0_m,u1_m,grad_matrix[id_thread],nbt,n0,tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWuWv[id_thread],m,freq_matrix[id_thread],nbk,v0_m,v1_m,tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread])

        for k=1:nbk

            tab_Mpq_c_bare_par_sampling_rot!(tab_Mpq_c_bare_par[id_thread],tab_W_temp[id_thread],tabWu_kn_temp[id_thread],tabWuJ_kn_temp[id_thread],tabWv_kl_temp[id_thread],tabWvJ_kl_temp[id_thread],tab_omega,nbomega,k,m,Ftot,Omegau,Omegav,Omegaz,jacu,jacv,nbn,nbp,n0)

        end
    end

    tab_Mpq_a_bare = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b_bare = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c_bare = zeros(ComplexF64, nbp, nbp, nbomega)

    pref  = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz
    prefc = 2.0*(2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv

    Threads.@threads for index=1:nbp*nbp 

        tab_Mpq_bare_sampling!(index,tab_Mpq_a_bare,tab_Mpq_a_bare_par,nbp,nbomega,pref)
        tab_Mpq_bare_sampling!(index,tab_Mpq_b_bare,tab_Mpq_b_bare_par,nbp,nbomega,pref)
        tab_Mpq_bare_sampling!(index,tab_Mpq_c_bare,tab_Mpq_c_bare_par,nbp,nbomega,prefc)

    end

    return tab_Mpq_a_bare, tab_Mpq_b_bare, tab_Mpq_c_bare, tab_omega
end



function ResponseMatrix_m_GS_sampling_rot_split_separate(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, 
    nbRe::Int64, nbIm::Int64,nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, 
    nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbomega = nbRe * nbIm 

    tab_Mpq_a_bare, tab_Mpq_b_bare, tab_Mpq_c_bare, tab_omega = ResponseMatrix_m_bare_sampling_rot_split_separate(m,re_omega_min,re_omega_max,im_omega_min,im_omega_max,nbRe,nbIm,nbJu,nbJv,nbLz,nbt,epsLz,iJmin,iJmax,iJmin2d,iJmax2d)

    tab_Mpq_a = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c = zeros(ComplexF64, nbp, nbp, nbomega)


    Threads.@threads for index=1:nbp*nbp 

        compute_tab_Mpq_sampling_rot!(m,n0,index,nbn,nbp,tab_Mpq_a,tab_Mpq_b,tab_Mpq_c,tab_Mpq_a_bare,tab_Mpq_b_bare,tab_Mpq_c_bare,nbomega)

    end

    return tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega

end