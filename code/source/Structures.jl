# https://www.matecdev.com/posts/julia-structs.html#mutable-structs
# using Setfield

# using SphericalHarmonics

mutable struct OrbitalElements

    E::Float64
    I3::Float64 
    Lz::Float64

    # Ju::Float64 
    # Jv::Float64

    u0::Float64 
    u1::Float64 

    v0::Float64 
    v1::Float64

    tab_Flmn::Vector{Float64} # reverse sampling
    # tab_uEff::Vector{Float64} # reverse sampling
    # tab_u::Vector{Float64} # reverse sampling
    # tab_xi::Vector{Float64} # reverse sampling
    tab_sinhu::Vector{Float64} # reverse sampling
    # tab_sinuEff::Vector{Float64} # reverse sampling
    tab_tpu::Vector{Float64} # reverse sampling

    tab_alpha_E::Vector{Float64}
    tab_alpha_I3::Vector{Float64}
    tab_alpha_Lz::Vector{Float64}




    tab_ylm::Vector{Float64} # reverse sampling
    # tab_vEff::Vector{Float64} # reverse sampling
    # tab_v::Vector{Float64} # reverse sampling
    tab_sinv::Vector{Float64} # reverse sampling
    # tab_sinvEff::Vector{Float64} # reverse sampling
    tab_tpv::Vector{Float64} # reverse sampling

    tab_beta_E::Vector{Float64}
    tab_beta_I3::Vector{Float64}
    tab_beta_Lz::Vector{Float64}

    # container_Y # https://jishnub.github.io/SphericalHarmonics.jl/stable/#SphericalHarmonics.cache-Tuple
    

    # alphak_E_2 = Delta^2 *sinhu^2/tpu
    #     alphak_I3_2 = -Delta^2 /tpu
    #     alphak_Lz_2 =



    # matrix_freq::Matrix{Float64}
    # matrix_grad::Matrix{Float64}


end

function OrbitalElements_init(nbt::Int64)

    # tab_uEff = zeros(Float64,2*nbt+1)
    # tab_u = zeros(Float64,2*nbt+1)
    # tab_xi = zeros(Float64,2*nbt+1)
    tab_sinhu = zeros(Float64,2*nbt+1)
    # tab_sinuEff = zeros(Float64,2*nbt+1)
    tab_tpu = zeros(Float64,2*nbt+1)
    tab_Flmn = zeros(Float64,2*nbt+1)

    tab_alpha_E = zeros(Float64,2*nbt+1)
    tab_alpha_I3 = zeros(Float64,2*nbt+1)
    tab_alpha_Lz = zeros(Float64,2*nbt+1)

    tab_ylm = zeros(Float64,2*nbt+1)
    # tab_vEff = zeros(Float64,2*nbt+1)
    # tab_v = zeros(Float64,2*nbt+1)
    tab_sinv = zeros(Float64,2*nbt+1)
    # tab_sinvEff = zeros(Float64,2*nbt+1)
    tab_tpv = zeros(Float64,2*nbt+1)

    tab_beta_E = zeros(Float64,2*nbt+1)
    tab_beta_I3 = zeros(Float64,2*nbt+1)
    tab_beta_Lz = zeros(Float64,2*nbt+1)

    # container_Y = SphericalHarmonics.cache(lmax) 

    elements = OrbitalElements(0.0,0.0,0.0,0.0,0.0,0.0,0.0,tab_Flmn,tab_sinhu,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_sinv,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)

    return elements

end


function OrbitalElements_update!(elem::OrbitalElements, E::Float64, I3::Float64, Lz::Float64, l::Int64, m::Int64, n::Int64, nbt::Int64)

    # Ju = _Ju(E,Lz,I3)
    # Jv = _Jv(E,Lz,I3)

    u0, u1 = find_bounds_u(E,Lz,I3)
    v0, v1 = find_bounds_v(E,Lz,I3)

    elem.E = E 
    elem.Lz = Lz 
    elem.I3 = I3
    # elem.Ju = Ju 
    # elem.Jv = Jv



    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

    # println("a0")
    
    for iv=0:2*nbt # reverse sampling
        # println("b0")
        vEff = pi/2 - pi/(2*nbt) * iv
        sinvEff = sin(vEff)
        
        v = tv*sinvEff + sv
        sinv = sin(v)

        # println("b1")

        # ylm1 = computePlmcostheta!(elem.container_Y, v, l)[(l,m)]
        
        # println(ylm)
        ylm = SphericalHarmonics.associatedLegendre(v,l,m)/sqrt(2.0)
        # ylm = plm/sqrt(2.0)
        # println((ylm1,ylm))

        # println("b2")
        # ylm = elem.container_Y[(l,m)]

        # println("b3")

        # ylm = Ylm(l,m,v,0.0)
        tpv =  _tpv(vEff,E,Lz,I3,v0,v1)

        # println(ylm)
        elem.tab_ylm[iv+1] = real(ylm)
        # elem.tab_vEff[iv+1] = vEff
        # elem.tab_v[iv+1] = v
        elem.tab_sinv[iv+1] = sinv
        # elem.tab_sinvEff[iv+1] = sinvEff

        elem.tab_tpv[iv+1] = tpv

        elem.tab_beta_E[iv+1] = Delta^2 *sinv^2/tpv
        elem.tab_beta_I3[iv+1] = Delta^2 /tpv
        elem.tab_beta_Lz[iv+1] = - Lz/(tpv*sinv^2)

 
    end

    # println("a1")

    for iu=0:2*nbt # reverse sampling
        # println("c0")
        uEff = pi/2 - pi/(2*nbt) * iu
        sinuEff = sin(uEff)
        u = tu*sinuEff + su
        sinhu = sinh(u)
        xi = cosh(u)
        Flmn = itp_Flmn_xi_exact(l,m,n,xi)
        
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)

        # elem.tab_uEff[iu+1] = uEff
        # elem.tab_u[iu+1] = u 
        # elem.tab_xi[iu+1] = xi
        elem.tab_sinhu[iu+1] = sinhu
        # elem.tab_sinuEff[iu+1] = sinuEff

        elem.tab_tpu[iu+1] = tpu

       
        
        elem.tab_Flmn[iu+1] = Flmn

        elem.tab_alpha_E[iu+1] = Delta^2 *sinhu^2/tpu
        elem.tab_alpha_I3[iu+1] = -Delta^2 /tpu
        elem.tab_alpha_Lz[iu+1] = -Lz/(tpu*sinhu^2)

 
    end

    # println("a2")
    # matrix_freq, matrix_grad = frequency_matrix(E,Lz,I3)

    # elements = OrbitalElements(E,I3,Lz,Ju,Jv,u0,u1,v0,v1,tab_Flmn,tab_uEff,tab_u,tab_xi,tab_sinhu,tab_sinuEff,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_vEff,tab_v,tab_sinv,tab_sinvEff,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)
    # elements = OrbitalElements(E,I3,Lz,u0,u1,v0,v1,tab_Flmn,tab_uEff,tab_u,tab_xi,tab_sinhu,tab_sinuEff,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_vEff,tab_v,tab_sinv,tab_sinvEff,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)

    # return elements

end