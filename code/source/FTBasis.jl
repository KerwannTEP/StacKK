# include("RadialBasisExact.jl")
# # include("RadialBasis.jl)
# include("SpheroidalHarmonics.jl")
# include("Args.jl")
# include("Constants.jl")

using LinearAlgebra

using Plots

# Backward integration
# USE RK4
# Tqke inspiration from https://github.com/JuliaStellarDynamics/OrbitalElements.jl/blob/main/src/Utils/Integrators.jl
function tab_alphak(k1::Int64, k2::Int64, k3::Int64, Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=100)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    u0, u1 = find_bounds_u(E,Lz,I3)
    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    # println((u0,u1))

    trans_freq = transpose(frequency_matrix(E,Lz,I3))
    trans_bk = [k1 k2 k3]

    # display(trans_freq)

    sum = [0.0, 0.0, 0.0]

    tab_alphak = zeros(Float64, 2, nbu) # (u, alphak[u] )

    uEff = pi/2
    du = pi/nbu

    # Test 

    tabint = zeros(Float64, 2, nbu)
    


    for i=1:nbu 
        # uEff = pi/2 - pi/nbu*(i-0.5)

        # Step 1
        u = tu*sin(uEff) + su
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)
        
        sum[1] += du/6.0 * Delta^2 *sinh(u)^2/tpu
        sum[2] += -du/6.0 * Delta^2 /tpu
        sum[3] += -du/6.0 * Lz/(tpu*sinh(u)^2)

        # Step 2 and 3
        uEff -= du/2.0
        u = tu*sin(uEff) + su
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)

        tabint[1,i],  tabint[2,i] = uEff+pi/2, Lz/(tpu*sinh(u)^2)

        sum[1] += 4.0*du/6.0 * Delta^2 *sinh(u)^2/tpu
        sum[2] += -4.0*du/6.0 * Delta^2 /tpu
        sum[3] += -4.0*du/6.0 * Lz/(tpu*sinh(u)^2)

         # Step 4
         uEff -= du/2.0
         u = tu*sin(uEff) + su
         tpu = _tpu(uEff,E,Lz,I3,u0,u1)
 
         sum[1] += du/6.0 * Delta^2 *sinh(u)^2/tpu
         sum[2] += -du/6.0 * Delta^2 /tpu
         sum[3] += -du/6.0 * Lz/(tpu*sinh(u)^2)


        # println((sum[1],sum[2],sum[3]))


 
        # display(trans_freq*sum)
        # display(bk)

        akr = trans_bk*trans_freq*sum
        
        tab_alphak[1,i], tab_alphak[2,i] = uEff, k1*pi-akr[1]
    end

    # pt = plot(tabint[1,:],  tabint[2,:], xaxis=:log10)
    # savefig(pt,"test_integrand.png")

    return tab_alphak

end

# Backward integration
# USE RK4
function tab_betak(k1::Int64, k2::Int64, k3::Int64, Ju::Float64, Jv::Float64, Lz::Float64, nbv::Int64=100)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    v0, v1 = find_bounds_v(E,Lz,I3)
    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

 

    trans_freq = transpose(frequency_matrix(E,Lz,I3))
    trans_bk = [k1 k2 k3]

    sum = [0.0, 0.0, 0.0]

    tab_betak = zeros(Float64, 2, nbv) # (v, betak[v] )


    for i=1:nbv 
        vEff = pi/2 - pi/nbv*(i-0.5)
        v = tv*sin(vEff) + sv
        tpv = _tpv(vEff,E,Lz,I3,v0,v1)
        

        sum[1] += pi/nbv * Delta^2 *sin(v)^2/tpv
        sum[2] += pi/nbv * Delta^2 /tpv
        sum[3] += -pi/nbv * Lz/(tpv*sin(v)^2)

 
        # display(trans_freq*sum)
        # display(bk)

        bkr = trans_bk*trans_freq*sum
        
        tab_betak[1,i], tab_betak[2,i] = vEff, k2*pi-bkr[1]
    end

    return tab_betak

end



# Use exact formalism here
# https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
# https://github.com/JuliaStellarDynamics/LinearResponse.jl/blob/main/src/WMat.jl
# function _Wklmn(k1::Int64, k2::Int64, k3::Int64, l::Int64, m::Int64, n::Int64, Ju::Float64, Jv::Float64, Lz::Float64, nbt::Int64=100)
function _Wklmn(k1::Int64, k2::Int64, k3::Int64, elem::OrbitalElements, trans_freq, matrix_grad::Matrix{Float64}, nbt::Int64=100)

    # Compute alpha_k, beta_k on the fly 
    # Backward integration


    Jw = abs(matrix_grad[1,1]*matrix_grad[2,2]-matrix_grad[2,1]*matrix_grad[1,2])

  
    Wu1 = 0.0 # With Jwu
    Wu2 = 0.0
    Wv1 = 0.0
    Wv2 = 0.0 # With Jwv

    # Backward integration 
    # RK4

    du = pi/nbt 
    alphak = k1*pi

    iu = 1 # index to access elem.tab_ylm
    sinhu = elem.tab_sinhu[iu]
    sinuEff = elem.tab_sinuEff[iu]
    xi = elem.tab_xi[iu]
    u = elem.tab_u[iu]
    uEff = elem.tab_uEff[iu]

    x = (xi-1)/(xi+1)
    tpu = elem.tab_tpu[iu]
    Jwu = Delta^4/Jw * sinhu^2

    Flmn = elem.tab_Flmn[iu]

    # pref_mat = trans_bk*trans_freq
    pref_mat_1 = k1*trans_freq[1,1] + k2*trans_freq[2,1] + k3*trans_freq[3,1] 
    pref_mat_2 = k1*trans_freq[1,2] + k2*trans_freq[2,2] + k3*trans_freq[3,2] 
    pref_mat_3 = k1*trans_freq[1,3] + k2*trans_freq[2,3] + k3*trans_freq[3,3] 

    for i=1:nbt 

        # Step 1

        pref_Wu2_1 = du/6.0 * 1.0/pi * Flmn*exp(-1im*alphak)/tpu
        pref_Wu1_1 = Jwu * pref_Wu2_1


        Wu1 += pref_Wu1_1
        Wu2 += pref_Wu2_1

        # Update alphak 
        alphak_E_1 = elem.tab_alpha_E[iu] 
        alphak_I3_1 = elem.tab_alpha_I3[iu] 
        alphak_Lz_1 = elem.tab_alpha_Lz[iu] 
        dalphak_1 = du*(pref_mat_1*alphak_E_1 + pref_mat_2*alphak_I3_1 + pref_mat_3*alphak_Lz_1)

        # Step 2

        iu += 1 # index to access elem.tab_ylm
        sinhu = elem.tab_sinhu[iu]
        sinuEff = elem.tab_sinuEff[iu]
        xi = elem.tab_xi[iu]
        u = elem.tab_u[iu]
        uEff = elem.tab_uEff[iu]

        x = (xi-1)/(xi+1)
        tpu = elem.tab_tpu[iu]
        Jwu = Delta^4/Jw * sinhu^2

        Flmn = elem.tab_Flmn[iu]

        pref_Wu2_2 = du/3.0 * 1/pi * Flmn*exp(-1im*(alphak-0.5*dalphak_1))/tpu
        pref_Wu1_2 =  Jwu*pref_Wu2_2#

        Wu1 += pref_Wu1_2
        Wu2 += pref_Wu2_2

        # Update alphak 
        alphak_E_2 = elem.tab_alpha_E[iu] 
        alphak_I3_2 = elem.tab_alpha_I3[iu]
        alphak_Lz_2 = elem.tab_alpha_Lz[iu]

        dalphak_2 = du*(pref_mat_1*alphak_E_2 + pref_mat_2*alphak_I3_2 + pref_mat_3*alphak_Lz_2)


        # Step 3

        pref_Wu2_3 = du/3.0 * 1/pi * Flmn*exp(-1im*(alphak-0.5*dalphak_2))/tpu
        pref_Wu1_3 = Jwu*pref_Wu2_3#
        
        Wu1 += pref_Wu1_3
        Wu2 += pref_Wu2_3

        # Update alphak 
        dalphak_3 = dalphak_2

        # Step 4

        iu += 1 # index to access elem.tab_ylm
        sinhu = elem.tab_sinhu[iu]
        sinuEff = elem.tab_sinuEff[iu]
        xi = elem.tab_xi[iu]
        u = elem.tab_u[iu]
        uEff = elem.tab_uEff[iu]

        x = (xi-1)/(xi+1)
        tpu = elem.tab_tpu[iu]
        Jwu = Delta^4/Jw * sinhu^2

        Flmn = elem.tab_Flmn[iu]

        pref_Wu2_4 = du/6.0 * 1/pi * Flmn*exp(-1im*(alphak-dalphak_3))/tpu
        pref_Wu1_4 = Jwu*pref_Wu2_4

        Wu1 += pref_Wu1_4
        Wu2 += pref_Wu2_4

        # Update alphak 
        alphak_E_4 = elem.tab_alpha_E[iu]
        alphak_I3_4 = elem.tab_alpha_I3[iu]
        alphak_Lz_4 = elem.tab_alpha_Lz[iu] 

        dalphak_4 = du*(pref_mat_1*alphak_E_4 + pref_mat_2*alphak_I3_4 + pref_mat_3*alphak_Lz_4)

        # Update alphak for next iteration
        alphak -= (dalphak_1 + 2.0*dalphak_2 + 2.0*dalphak_3 + dalphak_4)/6.0
         

    end

    dv = pi/nbt 
    betak = k2*pi

    iv = 1 # index to access elem.tab_ylm
    sinv = elem.tab_sinv[iv]
    sinvEff = elem.tab_sinvEff[iv]
    v = elem.tab_v[iv]
    vEff = elem.tab_vEff[iv]
   
    ylm = elem.tab_ylm[iv]
    tpv = elem.tab_tpv[iv]
    Jwv = Delta^4/Jw * sinv^2

    


    for i=1:nbt 

        # Step 1

        pref_Wv1_1 = dv/6.0 * 1/pi * ylm*exp(-1im*betak)/tpv
        pref_Wv2_1 = Jwv*pref_Wv1_1#

        Wv1 += pref_Wv1_1
        Wv2 += pref_Wv2_1

        # Update betak 

        betak_E_1 =  elem.tab_beta_E[iv] 
        betak_I3_1 = elem.tab_beta_I3[iv] 
        betak_Lz_1 = elem.tab_beta_Lz[iv] 

        dbetak_1 = dv*(pref_mat_1*betak_E_1 + pref_mat_2*betak_I3_1 + pref_mat_3*betak_Lz_1)

        # Step 2

        iv += 1 # update
        sinv = elem.tab_sinv[iv]
        sinvEff = elem.tab_sinvEff[iv]
        v = elem.tab_v[iv]
        vEff = elem.tab_vEff[iv]

        tpv = elem.tab_tpv[iv]
        Jwv = Delta^4/Jw * sinv^2
        ylm = elem.tab_ylm[iv]

        pref_Wv1_2 = dv/3.0 * 1/pi * ylm*exp(-1im*(betak-0.5*dbetak_1))/tpv
        pref_Wv2_2 = Jwv*pref_Wv1_2

        Wv1 += pref_Wv1_2
        Wv2 += pref_Wv2_2

        # Update betak 
        betak_E_2 = elem.tab_beta_E[iv]
        betak_I3_2 =  elem.tab_beta_I3[iv]
        betak_Lz_2 = elem.tab_beta_Lz[iv] 

        dbetak_2 = dv*(pref_mat_1*betak_E_2 + pref_mat_2*betak_I3_2 + pref_mat_3*betak_Lz_2)


        # Step 3

        pref_Wv1_3 = dv/3.0 * 1/pi * ylm*exp(-1im*(betak-0.5*dbetak_2))/tpv
        pref_Wv2_3 = Jwv*pref_Wv1_3
        
        Wv1 += pref_Wv1_3
        Wv2 += pref_Wv2_3

        # Update betak 
        dbetak_3 = dbetak_2

        # Step 4
        iv += 1 # update 
        sinv = elem.tab_sinv[iv]
        sinvEff = elem.tab_sinvEff[iv]
        v = elem.tab_v[iv]
        vEff = elem.tab_vEff[iv]

        tpv = elem.tab_tpv[iv]
        Jwv = Delta^4/Jw * sinv^2
        ylm = elem.tab_ylm[iv]


        pref_Wv1_4 = dv/6.0 * 1/pi * ylm*exp(-1im*(betak-dbetak_3))/tpv
        pref_Wv2_4 = Jwv*pref_Wv1_4

        Wv1 += pref_Wv1_4
        Wv2 += pref_Wv2_4



        # Update betak 
        betak_E_4 = elem.tab_beta_E[iv] 
        betak_I3_4 =  elem.tab_beta_I3[iv] 
        betak_Lz_4 = elem.tab_beta_Lz[iv] 

        dbetak_4 = dv*(pref_mat_1*betak_E_4 + pref_mat_2*betak_I3_4 + pref_mat_3*betak_Lz_4)

        # Update betak for next iteration
        betak -= (dbetak_1 + 2.0*dbetak_2 + 2.0*dbetak_3 + dbetak_4)/6.0
         
    end

    return Wu1*Wv1 + Wu2*Wv2 

end



function ResponseMatrixElement(omega::ComplexF64, lp::Int64, mp::Int64, np::Int64, lq::Int64, mq::Int64, nq::Int64, nbJ::Int64=100, nbt::Int64=100)

    if (mp != mq)
        return 0.0

    else 

        k3 = mp # also equal to mq
        sum = 0.0

        elem_p = OrbitalElements_init(nbt)
        elem_q = OrbitalElements_init(nbt)

        freq_matrix = zeros(Float64, 3, 3)
        grad_matrix = zeros(Float64, 3, 3)


       for iu=1:nbJ 
            Ju = Jumax/nbJ*(iu-0.5)
            for iv=1:nbJ 
                Jv = Jvmax/nbJ*(iv-0.5)
                for iz=1:nbJ 
                    Lz = -Lzmax + 2.0*Lzmax/nbJ*(iz-0.5)

                    # println((Ju,Jv,Lz))

                    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

                    dFdE = _dFdE(E,Lz)
                    dFdLz = _dFdLz(E,Lz)

                 
                    # precompute table of Ylm(l,m,v,0,0) here

                    OrbitalElements_update!(elem_p,E,I3,Lz,lp,mp,np,nbt)
                    OrbitalElements_update!(elem_q,E,I3,Lz,lq,mq,nq,nbt)


                    fill_grad_frequency_matrix!(grad_matrix,freq_matrix,E,Lz,I3)


                     trans_freq = transpose(freq_matrix)
                    Omegau, Omegav, Omegaz = freq_matrix[1,1], freq_matrix[1,2], freq_matrix[1,3] # TODO 

                    dFdJu = dFdE * Omegau # dF/dJu = dF/dE dE/dJu + dF/dI3 dI3/dJu + dF/dLz dLz/dJu
                    dFdJv = dFdE * Omegav # dF/dJv = dF/dE dE/dJv + dF/dI3 dI3/dJv + dF/dLz dLz/dJv
                    dFdLz = dFdE * Omegaz + dFdLz # dF/dLz = dF/dE dE/dLz + dF/dI3 dI3/dLz + dF/dLz dLz/dLz

                    
                    # println("====")

                    for k1=-kmax:kmax 
                        # println(k1)
                        for k2=-kmax:kmax
                            # println("k2=",k2)


                            # Fill directly all matrix elements here, in this loop 
                            # This is because we don't want to repeat some of the computation before
                            # Only Wkp, Wkq depends on p, q 
                            # Threading on this loop 
                            # Fixed m
                            # If we have (p, q) then we have (q, p) immediately
                            # So we can cut in half the computation
                            # Loop on l begins at lp = m and lq = m
                            # Loop on n begins at 0 for m=0 and at 1 for m != 0
                            # Are dodes for -m are related to the modes for m ? (look into this)

                            # println("c0")
                            Wkp = _Wklmn(k1,k2,k3,elem_p,trans_freq, grad_matrix,nbt)
                            # println("c1")
                             Wkq = _Wklmn(k1,k2,k3,elem_q,trans_freq, grad_matrix,nbt)

                            # println((Wkp,Wkq))
                            # println("c2")

                            kdotdF = k1*dFdJu + k2*dFdJv + k3*dFdLz

                            
                            kdotOmega = k1*Omegau + k2*Omegav + k3*Omegaz 

                            # println((omega-kdotOmega))

                            sum += kdotdF/(omega-kdotOmega) * conj(Wkp) * Wkq 
                            # println((Ju,Jv,sum))

                            # println("-=-=-=-")

                        end
                        
                    end
                    # println("----")

                    # println((Omegau, Omegav, Omegaz))
                end
                # println((Ju,Jv,sum))

            end
        end

        sum *= (2*pi)^3 * sqrt(4*pi*G)/Delta * Jumax/nbJ * Jvmax/nbJ * 2.0*Lzmax/nbJ

        return sum 
    end
end

function ResponseMatrix_m(m::Int64, omega::ComplexF64, nbJ::Int64=100, nbt::Int64=100)

    # Compute ResponseMatrix for harmonics m
    # Indexes are p=(l,n) at fixed m.
    # n between n0 and nmax: n0=0 for m=0 // n0=1 for m=1

    @assert (abs(m) <= mmax) "m must be less that m_max"

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbp = (lmax - abs(m) + 1) * (nmax - n0 + 1)
    

    tab_pq = zeros(Int64, nbp*nbp, 2)
    tab_lpnplqnq = zeros(Int64, nbp*nbp, 4)

    index = 1
    p = 1
    for lp=abs(m):lmax 
        for np=n0:nmax 
            q = 1
            for lq=abs(m):lmax 
                for nq=n0:nmax 

                    tab_lpnplqnq[index,1], tab_lpnplqnq[index,2], tab_lpnplqnq[index,3], tab_lpnplqnq[index,4] = lp, np, lq, nq

                    tab_pq[index,1], tab_pq[index,2] = p, q

                    index += 1
                    q += 1
                end
            end

            p += 1
        end
    end

    tab_Mpq_par = zeros(ComplexF64, nbp*nbp)

    tab_Mpq = zeros(ComplexF64, nbp, nbp)


    progress = Threads.Atomic{Int64}(0)
    nbgrid = nbp*nbp

    Threads.@threads for index=1:nbgrid
    # for index=1:nbgrid


        lp, np, lq, nq = tab_lpnplqnq[index,1], tab_lpnplqnq[index,2], tab_lpnplqnq[index,3], tab_lpnplqnq[index,4]
        Mpq = ResponseMatrixElement(omega,lp,m,np,lq,m,nq,nbJ,nbt)

        tab_Mpq_par[index] = Mpq


        Threads.atomic_add!(progress,1)

        println("Progress = ", progress[] ," / ",nbp*nbp) #, " (p,q)=",(p,q)," (lp,np,lq,nq)=",(lp,np,lq,nq), " Mpq=", Mpq)


    end

    for index=1:nbp*nbp 

        p, q = tab_pq[index,1], tab_pq[index,2]
        tab_Mpq[p,q] = tab_Mpq_par[index]

    end

    tab_ln = zeros(Int64, nbp, 2) # (l,n)

    p = 1

    for l=abs(m):lmax 
        for n=n0:nmax 

            tab_ln[p,1], tab_ln[p,2] = l, n 
            p += 1
        end
    end

    return tab_Mpq, tab_ln


end

function det_Dielectric(tabMpq, size)

    id = 1.0* Matrix(I, size, size)

    return det(id - tabMpq)

end

