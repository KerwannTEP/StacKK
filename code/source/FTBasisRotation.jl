
############################################################################################
# Sample the response matrix over a grid of complex frequencies
############################################################################################

# compute Mpq_a, Mpq_b, Mpq_c
function ResponseMatrix_m_sampling_rot(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default, epsLz::Float64=4.0*10^(-5))

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



    tab_Mpq_a_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]
    tab_Mpq_b_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()]

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()]


    # 3D integration for integralds a and b
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

                        tab_Mpq_a_par[index,iomega,id_thread] += integrand * jacu * jacv * jacz
                        tab_Mpq_b_par[index,iomega,id_thread] += integrand * sign(Lz) * jacu * jacv * jacz
                    end

                end

            end

        end

        # println("-----")


    end



    # 2D integration for c
    # Evaluate at Lz = +epsLz and Lz = -epsLz


    nbgrid2d = nbJ^2
    tab_Mpq_c_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]


    Threads.@threads for iJ=1:nbgrid2d

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

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_p,I3_p,Lz_p,u0_p,u1_p,v0_p,v1_p,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, positive part
                    end

                end

            end

        end




        # Lz = 0-

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_m,Lz_m,I3_m)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_m,Lz_m)

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_m,I3_m,Lz_m,u0_m,u1_m,v0_m,v1_m,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, negative part
                    end

                end

            end

        end




    end









    tab_Mpq_a = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c = zeros(ComplexF64, nbp, nbp, nbomega)


    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz


    Threads.@threads for index=1:nbp*nbp

        # index - 1 = (q-1) + (p-1)*nbp

        p = floor(Int64, (index-1)/nbp) + 1
        q  = index  - (p-1)*nbp

        for iomega=1:nbomega


            # tab_Mpq[p,q,iomega] = tab_Mpq_par_re[index,iomega][] + 1im*tab_Mpq_par_im[index,iomega][]

            for it=1:Threads.nthreads()

                # tab_Mpq[p,q,iomega] += tab_Mpq_par_re_par[index,iomega,it] + 1im*tab_Mpq_par_im_par[index,iomega,it]
                tab_Mpq_a[p,q,iomega] +=tab_Mpq_a_par[index,iomega,it]
                tab_Mpq_b[p,q,iomega] +=tab_Mpq_b_par[index,iomega,it]
                tab_Mpq_c[p,q,iomega] +=tab_Mpq_c_par[index,iomega,it]


            end
            tab_Mpq_a[p,q,iomega] *= pref
            tab_Mpq_b[p,q,iomega] *= pref
            tab_Mpq_c[p,q,iomega] *= 2.0*(2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv
        end
    end

    # pt=plot([tab_integrand[1:nbJ_default,1] tab_integrand[1:nbJ_default,2]])
    # savefig(pt, "integrand.png")

    return tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega
end



############################################################################################
# Sample the response matrix over a grid of complex frequencies
# Split the action integral to distribute the load over multiple nodes
############################################################################################

function ResponseMatrix_m_sampling_rot_split(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, nbJ::Int64=nbJ_default, nbt::Int64=nbt_default, epsLz::Float64=4.0*10^(-5), iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

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

    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJ
    dtJv = pi/2.0*(1.0-epsJ)/nbJ
    dtLz = pi*(1.0-epsJ)/nbJ

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



    tab_Mpq_a_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]
    tab_Mpq_b_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()]

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()]


    # Maybe use different nbJu, nbJv, nbLz
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
        tLz = -pi/2.0+epsJ + dtLz*(iz-0.5)

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

                        tab_Mpq_a_par[index,iomega,id_thread] += integrand * jacu * jacv * jacz
                        tab_Mpq_b_par[index,iomega,id_thread] += integrand * sign(Lz) * jacu * jacv * jacz
                    end

                end

            end

        end

        # println("-----")


    end

    # 2D integration for c
    # Evaluate at Lz = +epsLz and Lz = -epsLz


    nbgrid2d = nbJ^2
    tab_Mpq_c_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]




    # Use iJmin2d and iJmax2d

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

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_p,I3_p,Lz_p,u0_p,u1_p,v0_p,v1_p,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, positive part
                    end

                end

            end

        end




        # Lz = 0-

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_m,Lz_m,I3_m)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_m,Lz_m)

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_m,I3_m,Lz_m,u0_m,u1_m,v0_m,v1_m,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, negative part
                    end

                end

            end

        end




    end







    tab_Mpq_a = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c = zeros(ComplexF64, nbp, nbp, nbomega)


    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz


    Threads.@threads for index=1:nbp*nbp

        # index - 1 = (q-1) + (p-1)*nbp

        p = floor(Int64, (index-1)/nbp) + 1
        q  = index  - (p-1)*nbp

        for iomega=1:nbomega


            # tab_Mpq[p,q,iomega] = tab_Mpq_par_re[index,iomega][] + 1im*tab_Mpq_par_im[index,iomega][]

            for it=1:Threads.nthreads()

                # tab_Mpq[p,q,iomega] += tab_Mpq_par_re_par[index,iomega,it] + 1im*tab_Mpq_par_im_par[index,iomega,it]
                tab_Mpq_a[p,q,iomega] +=tab_Mpq_a_par[index,iomega,it]
                tab_Mpq_b[p,q,iomega] +=tab_Mpq_b_par[index,iomega,it]
                tab_Mpq_c[p,q,iomega] +=tab_Mpq_c_par[index,iomega,it]


            end
            tab_Mpq_a[p,q,iomega] *= pref
            tab_Mpq_b[p,q,iomega] *= pref
            tab_Mpq_c[p,q,iomega] *= 2.0*(2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv

        end
    end

    # pt=plot([tab_integrand[1:nbJ_default,1] tab_integrand[1:nbJ_default,2]])
    # savefig(pt, "integrand.png")

    return tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega
end


function det_Dielectric_rot(A, B, C, alpha, size)

    id = 1.0* Matrix(I, size, size)
    M = A + alpha*(B+C)

    return det(id - M)

end


############################################################################################
# Sample the response matrix over a grid of complex frequencies
# Separate Ju,Jv action space sampling from Lz action space sampling
############################################################################################


function ResponseMatrix_m_sampling_rot_split_separate(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, 
                                nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), 
                                iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

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

    epsJ = 0.01 # Action cutoff

    dtJu = pi/2.0*(1.0-epsJ)/nbJu
    dtJv = pi/2.0*(1.0-epsJ)/nbJv
    dtLz = pi*(1.0-epsJ)/nbLz

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


    elem = [OrbitalElements_alt_half_init(nbt) for it=1:Threads.nthreads()]

    freq_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]
    grad_matrix = [zeros(Float64, 3, 3)  for it=1:Threads.nthreads()]



    tab_Mpq_a_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]
    tab_Mpq_b_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]

    tabWkp_temp = [zeros(Float64, nbp, nbk) for it=1:Threads.nthreads()]

    tabWuWv = [zeros(Float64, nbk, 4) for it=1:Threads.nthreads()]


    # Maybe use different nbJu, nbJv, nbLz
    Threads.@threads for iJ=iJmin:iJmax

        # iJ - 1 = iz-1 + nbz*(iv-1) + nbz*nbv*(iu-1)
        # nbz = nbJ
        # nbv = nbJ

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

                        tab_Mpq_a_par[index,iomega,id_thread] += integrand * jacu * jacv * jacz
                        tab_Mpq_b_par[index,iomega,id_thread] += integrand * sign(Lz) * jacu * jacv * jacz
                    end

                end

            end

        end

        # println("-----")


    end

    # 2D integration for c
    # Evaluate at Lz = +epsLz and Lz = -epsLz


    tab_Mpq_c_par = [0.0 + 0.0im for index=1:nbp*nbp, iomega=1:nbomega, it=1:Threads.nthreads()]




    # Use iJmin2d and iJmax2d

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

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_p,I3_p,Lz_p,u0_p,u1_p,v0_p,v1_p,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, positive part
                    end

                end

            end

        end




        # Lz = 0-

        fill_grad_frequency_matrix!(grad_matrix[id_thread],freq_matrix[id_thread],E_m,Lz_m,I3_m)
        Omegau, Omegav, Omegaz = freq_matrix[id_thread][1,1], freq_matrix[id_thread][1,2], freq_matrix[id_thread][1,3]

        Ftot = F(E_m,Lz_m)

        for p=1:nbp

            # nbn = nmax - n0 + 1
            # p = n-n0 + nbn*(l-abs(m)) + 1
            # p - 1 = n-n0 + nbn * l

            l = abs(m) + floor(Int64,(p-1)/nbn)
            n = n0 + p - 1 - nbn * (l - abs(m))



            OrbitalElements_alt_half_update!(elem[id_thread],l,m,n,E_m,I3_m,Lz_m,u0_m,u1_m,v0_m,v1_m,grad_matrix[id_thread],nbt)
            tabWkp_alt_half(p,l,tabWkp_temp[id_thread],tabWuWv[id_thread],m,elem[id_thread],freq_matrix[id_thread],nbk,nbt)


        end

        for k=1:nbk

            # nbkline = 2*kmax+1
            # k - 1 = (k2+kmax) + (k1+kmax)*(2*kmax+1)

            k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
            k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

            mF = m*Ftot

            kdotOmega = k1*Omegau + k2*Omegav + m*Omegaz

            for iomega=1:nbomega

                re_omega, im_omega = tab_omega[iomega,1], tab_omega[iomega,2]
                omega = re_omega + 1.0im*im_omega
                invden = mF/(omega-kdotOmega)

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

                        tab_Mpq_c_par[index,iomega,id_thread] += 0.5*integrand * jacu * jacv  # symmetric central value, negative part
                    end

                end

            end

        end




    end







    tab_Mpq_a = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_b = zeros(ComplexF64, nbp, nbp, nbomega)
    tab_Mpq_c = zeros(ComplexF64, nbp, nbp, nbomega)


    pref = (2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv*dtLz


    Threads.@threads for index=1:nbp*nbp

        # index - 1 = (q-1) + (p-1)*nbp

        p = floor(Int64, (index-1)/nbp) + 1
        q  = index  - (p-1)*nbp

        for iomega=1:nbomega


            # tab_Mpq[p,q,iomega] = tab_Mpq_par_re[index,iomega][] + 1im*tab_Mpq_par_im[index,iomega][]

            for it=1:Threads.nthreads()

                # tab_Mpq[p,q,iomega] += tab_Mpq_par_re_par[index,iomega,it] + 1im*tab_Mpq_par_im_par[index,iomega,it]
                tab_Mpq_a[p,q,iomega] +=tab_Mpq_a_par[index,iomega,it]
                tab_Mpq_b[p,q,iomega] +=tab_Mpq_b_par[index,iomega,it]
                tab_Mpq_c[p,q,iomega] +=tab_Mpq_c_par[index,iomega,it]


            end
            tab_Mpq_a[p,q,iomega] *= pref
            tab_Mpq_b[p,q,iomega] *= pref
            tab_Mpq_c[p,q,iomega] *= 2.0*(2*pi)^3 * 4*pi*G/Delta^2 *dtJu*dtJv

        end
    end

    # pt=plot([tab_integrand[1:nbJ_default,1] tab_integrand[1:nbJ_default,2]])
    # savefig(pt, "integrand.png")

    return tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega
end



############################################################################################
# Evaluate the response matrix at a complex frequency
# Separate Ju,Jv action space sampling from Lz action space sampling
############################################################################################


# function ResponseMatrix_m_sampling_rot_split_separate(m::Int64, re_omega_min::Float64, re_omega_max::Float64, im_omega_min::Float64, im_omega_max::Float64, nbRe::Int64, nbIm::Int64, 
#     nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3), 
#     iJmin::Int64=iJmin_default, iJmax::Int64=iJmax_default, iJmin2d::Int64=iJmin2d_default, iJmax2d::Int64=iJmax2d_default)

# @assert (abs(m) <= mmax) "m must be less that m_max"



function ResponseMatrix_m_rot_separate(m::Int64, omega::ComplexF64,  nbJu::Int64=nbJ_default, nbJv::Int64=nbJ_default, nbLz::Int64=80, 
                                                nbt::Int64=nbt_default, epsLz::Float64=1.0*10^(-3))

    tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, _ = ResponseMatrix_m_sampling_rot_split_separate(m, real(omega), real(omega), imag(omega), imag(omega), 1, 1, 
                                                                                                nbJu, nbJv, nbLz, nbt, epsLz, 1, nbJu*nbJv*nbLz, 1, nbJu*nbJv)

    return tab_Mpq_a[:,:,1], tab_Mpq_b[:,:,1], tab_Mpq_c[:,:,1]
end
