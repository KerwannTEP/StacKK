include("../source/Main.jl")




# Scalar product test (on x variable)
function scalar_productY(l::Int64, m::Int64, n::Int64, k::Int64, nbx::Int64=100)

    sum = 0.0

    for i=1:nbx
        x = 1.0/nbx*(i-0.5)
        Fx =  Flmn_x(l,m,n,x)
        Dx =  Dlmn_x(l,m,k,x)

        sum += 2/(1-x)^2*Fx*Dx
    end

    sum *= -(1.0/nbx)/Delta

    return sum
end

# Test Ylm(v,0.0)
function test_Ylm(l::Int64, m::Int64, v::Float64)

    ylm = Ylm(l,m,v,0.0)
    println("Ylm(v,0) [fct] = ",ylm)


    ylm_plm = SphericalHarmonics.associatedLegendre(v,l,m)/sqrt(2.0)
    println("Ylm(v,0) [Plm] = ",ylm_plm)

end







function test_fourier_psilmn(l::Int64, m::Int64, n::Int64, uEff::Float64, vEff::Float64, Ju::Float64, Jv::Float64, Lz::Float64, kmaxCV::Int64, nbt::Int64=nbt_default)


    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)
    u0, u1 = find_bounds_u(E,Lz,I3)
    v0, v1 = find_bounds_v(E,Lz,I3)

    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)
    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

    u = tu*sin(uEff) + su
    v = tv*sin(vEff) + sv

    # theory
    println("kmax = ",kmaxCV)

    xi = cosh(u)
    psith = Flmn(l,m,n,xi)*Ylm(l,m,v,0.0)
    println("psi (th)     = ",psith)

    # FT inverse
    # pu > 0 ; pv > 0

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)
    sump = [0.0 + 0.0im for it=1:Threads.nthreads()]
    sumprev = [0.0 + 0.0im for it=1:Threads.nthreads()]
    nbk = (2*kmax+1)
    nbk2 = nbk^2

    # Optimized

    n0 = 0
    if (m != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(m) + 1
    nbp = nbl * nbn
    nbk = (2*kmax+1)^2
    p = 1 + (l-abs(m))*nbn + (n-n0)
    # k = 1 + (k2+kmax) + (k1+kmax)*(2*kmax+1)

    elem = OrbitalElements_alt_half_init(nbt)
    freq_matrix = zeros(Float64, 3, 3)
    grad_matrix = zeros(Float64, 3, 3)
    tabWkp_temp = zeros(Float64, nbp, nbk)
    tabWuWv = zeros(Float64, nbk, 4)
    u0, u1 = find_bounds_u(E,Lz,I3)
    v0, v1 = find_bounds_v(E,Lz,I3)
    fill_grad_frequency_matrix!(grad_matrix,freq_matrix,E,Lz,I3)

    OrbitalElements_alt_half_update!(elem,l,m,n,E,I3,Lz,u0,u1,v0,v1,grad_matrix,nbt)
    tabWkp_alt_half(p,l,tabWkp_temp,tabWuWv,m,elem,freq_matrix,nbk,nbt)

    Threads.@threads for k=1:nbk2

        k1 =  -kmax + floor(Int64,(k-1)/(2*kmax+1))
        k2 = k - 1 - (k1+kmax)*(2*kmax+1) - kmax

        if ((abs(k1)<=kmaxCV) && (abs(k2)<=kmaxCV))

            it = Threads.threadid()

            alphak = compute_alphak(u,k1,k2,m,Ju,Jv,Lz)
            betak = compute_betak(v,k1,k2,m,Ju,Jv,Lz)

            if (mod(l+m+k2,2) ==0)

                Wkp = tabWkp_temp[p,k]

                sump[it] += Wkp*exp(1.0im*alphak)*exp(1.0im*betak)
                sumprev[it] += Wkp*exp(-1.0im*alphak)*exp(1.0im*betak)
            end
        end
        # end
    end

    sum = 0.0 + 0.0im
    sumrev = 0.0 + 0.0im
    for it=1:Threads.nthreads()
        sum += sump[it]
        sumrev += sumprev[it]
    end

    println("psi     (FT) = ",sum)
    println("psi Rev (FT) = ",sumrev)

end




# Compare psi(u,v) = (U[u]-V[v])/(sinh(u)^2+sin(v)^2)
function test_psi_Staeckel(u::Float64, v::Float64)

    psiStaeckel =  (_U(u)-_V(v))/(sinh(u)^2+sin(v)^2)
    lambda = lambda_from_u(u)
    nu = nu_from_v(v)
    psith = psi(lambda,nu)

    println("psi [theory]   = ",psith)
    println("psi [Staeckel] = ",psiStaeckel)

end






# U(r,r') = -G/|r-r'|
# (r-r')^2 = (x-x')^2 + (y-y')^2 + (z-z')^2
# x = R cos phi
# y = R sin phi
# R = Delta sinh u sin v
# z = Delta cosh u cos v
function test_interaction_potential(u,v,phi,up,vp,phip)

    xi = cosh(u)
    xip = cosh(up)

    R = Delta * sinh(u) * sin(v)
    z = Delta * cosh(u) * cos(v)
    Rp = Delta * sinh(up) * sin(vp)
    zp = Delta * cosh(up) * cos(vp)

    x = R * cos(phi)
    y = R * sin(phi)
    xp = Rp * cos(phip)
    yp = Rp * sin(phip)

    diff_r = sqrt((x-xp)^2 + (y-yp)^2 + (z-zp)^2)
    println("U(r,r') [th.] = ",-G/diff_r)

    sum = 0.0

    nblm = length(tab_lm)

    for lm=1:nblm
        l, m = tab_lm[lm]
        n0 = 0
        if (m != 0)
            n0 = 1
        end

        for n=n0:nmax

            # m >= 0
            psip_r = sqrt(4*pi*G)/Delta*Flmn(l,m,n,xi)*Ylm(l,m,v,phi)
            psip_rp = sqrt(4*pi*G)/Delta*Flmn(l,m,n,xip)*Ylm(l,m,vp,phip)

            sum += conj(psip_r)*psip_rp

            if (m != 0)

                psip_r = sqrt(4*pi*G)/Delta*Flmn(l,-m,n,xi)*Ylm(l,-m,v,phi)
                psip_rp = sqrt(4*pi*G)/Delta*Flmn(l,-m,n,xip)*Ylm(l,-m,vp,phip)

                sum += conj(psip_r)*psip_rp

            end
        end
    end

    sum *= -1.0

    println("U(r,r') [SCF] = ",sum)

end

# Compare rho(r) with DF integration
function test_DF_rho(u::Float64, v::Float64, nbE::Int64=100, nbv::Int64=100, nbK::Int64=100)

    lambda, nu = lambda_nu_from_u_v(u,v)
    R, z = R_z_from_lambda_nu(lambda,nu)
    psiRz = psi(lambda,nu)
    rhoRz = rho(lambda,nu)

    # Exact expression
    println("rho(u,v) = ",rhoRz)

    # Exact with R,z
    num = (a^2+c^2)*R^2 + 2*a^2*z^2 + 2*a^2*c^2+a^4+3*a^2*sqrt(a^2*c^2+c^2*R^2+a^2*z^2)
    den = (a^2*c^2+c^2*R^2+a^2*z^2)^(3/2)*(R^2+z^2+a^2+c^2+2*sqrt(a^2*c^2+c^2*R^2+a^2*z^2))^(3/2)
    println("rho(R,z) = ",M*c^2/(4*pi)*num/den)

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

    println("rho DF   = ",sum)
end



function test_dFdE(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    E_p = E + eps
    E_m = E - eps

    dFdE_num = (F(E_p,Lz,nbK) - F(E_m,Lz,nbK))/(2*eps)
    dFdE_th = _dFdE(E,Lz,nbK)

    println("dFdE Num = ",dFdE_num)
    println("dFdE Th  = ",dFdE_th)




end


function test_dFdLz(E::Float64, Lz::Float64, eps::Float64= 0.01, nbK::Int64=100)

    Lz_p = Lz + eps
    Lz_m = Lz - eps

    dFdLz_num = (F(E,Lz_p,nbK) - F(E,Lz_m,nbK))/(2*eps)
    dFdLz_th = _dFdLz(E,Lz,nbK)

    println("dFdLz Num = ",dFdLz_num)
    println("dFdLz Th  = ",dFdLz_th)


end



function test_dFdJ(Ju::Float64, Jv::Float64, Lz::Float64, dJ::Float64=0.001, nbK::Int64=100, nbu::Int64=100, nbv::Int64=100, eps::Float64=0.01)

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



    # Num

    # Ju
    E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju+dJ,Lz,Jv)
    E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju-dJ,Lz,Jv)

    Fp = F(E_p,Lz_p,nbK)
    Fm = F(E_m,Lz_m,nbK)

    dFdJu_num = (Fp-Fm)/(2.0*dJ)

     # Jv
     E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv+dJ)
     E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv-dJ)

     Fp = F(E_p,Lz_p,nbK)
     Fm = F(E_m,Lz_m,nbK)

     dFdJv_num = (Fp-Fm)/(2.0*dJ)

     # Lz
     E_p, Lz_p, I3_p = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz+dJ,Jv)
     E_m, Lz_m, I3_m = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz-dJ,Jv)

     Fp = F(E_p,Lz_p,nbK)
     Fm = F(E_m,Lz_m,nbK)

     dFdLz_num = (Fp-Fm)/(2.0*dJ)

     println("dFdJu (th ) = ",dFdJu)
     println("dFdJu (num) = ",dFdJu_num)
     println("------------------------")
     println("dFdJv (th ) = ",dFdJv)
     println("dFdJv (num) = ",dFdJv_num)
     println("------------------------")
     println("dFdLz (th ) = ",dFdLz)
     println("dFdLz (num) = ",dFdLz_num)

end


function frequency_test(Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=100, nbv::Int64=100, eps::Float64=0.001, dJ::Float64=0.001)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    dJudE, dJudI3, dJudLz = dJu(E,Lz,I3,nbu,eps)
    dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3,nbv)

   detJ = dJudE*dJvdI3-dJudI3*dJvdE

   Omega_u = dJvdI3/detJ
   Omega_v = -dJudI3/detJ
   Omega_z = (dJudI3*dJvdLz-dJudLz*dJvdI3)/detJ


   # num

   # Omegau = dEdJu

   E_p, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju+dJ,Lz,Jv)
   E_m, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju-dJ,Lz,Jv)

   Omega_u_num  = (E_p-E_m)/(2.0*dJ)

    # Omegav = dEdJv

    E_p, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv+dJ)
    E_m, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv-dJ)

    Omega_v_num  = (E_p-E_m)/(2.0*dJ)

     # Omegaz = dEdLz

     E_p, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz+dJ,Jv)
     E_m, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz-dJ,Jv)

     Omega_z_num  = (E_p-E_m)/(2.0*dJ)

     println("Omegau (th ) = ",Omega_u)
     println("Omegau (num) = ",Omega_u_num)
     println("------------------------")
     println("Omegav (th ) = ",Omega_v)
     println("Omegav (num) = ",Omega_v_num)
     println("------------------------")
     println("Omegaz (th ) = ",Omega_z)
     println("Omegaz (num) = ",Omega_z_num)

end


# integrate rho(xi,eta)
function integrate_rho(nbx::Int64=100, nbeta::Int64=100)

    sum = 0.0

    # int \rd \br rho(\br) = \Delta int \rd \xi \rd \eta \rd \phi c(\xi,\eta) \rho(\xi,\eta)
    # x = (xi-1)/(xi+1) in [0, 1]
    # xi = (1+x)/(1-x)
    # dxi = 2/(1-x)^2 dx

    # c(xi,eta) = Delta^2 xi^2 - Delta^2 eta^2

    for i=1:nbx
        x = 1.0/nbx*(i-0.5)
        xi = (1+x)/(1-x)
        lambda = c^2 + Delta^2*xi^2

        sumx = 0.0

        for j=1:nbeta
            eta = -1.0 + 2.0/nbeta*(j-0.5)
            nu = c^2 + Delta^2*eta^2
            v = acos(eta)

            rhor = rho(lambda,nu)
            cxieta = Delta^2*(xi^2 - eta^2)

            # println((rhor,cxieta,2.0/(1-x)^2))

            sumx += 2.0/nbeta * cxieta * rhor
        end

        sum += 1.0/nbx * 2.0/(1-x)^2 * sumx
    end

    sum *= 2.0*pi # integration over \phi yields 2*pi

    sum *= Delta

    println("Integration rho(r) = ",sum)
end
