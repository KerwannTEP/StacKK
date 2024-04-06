# https://www.matecdev.com/posts/julia-structs.html#mutable-structs
# using Setfield

# using SphericalHarmonics

mutable struct OrbitalElements

    E::Float64
    I3::Float64 
    Lz::Float64

    u0::Float64 
    u1::Float64 

    v0::Float64 
    v1::Float64

    tab_Flmn::Vector{Float64} # reverse sampling
    tab_xi::Vector{Float64} # reverse sampling
    tab_sinhu::Vector{Float64} # reverse sampling
    tab_tpu::Vector{Float64} # reverse sampling

    tab_alpha_E::Vector{Float64}
    tab_alpha_I3::Vector{Float64}
    tab_alpha_Lz::Vector{Float64}

    tab_ylm::Vector{Float64} # reverse sampling
    tab_v::Vector{Float64} # reverse sampling
    tab_sinv::Vector{Float64} # reverse sampling
    tab_tpv::Vector{Float64} # reverse sampling

    tab_beta_E::Vector{Float64}
    tab_beta_I3::Vector{Float64}
    tab_beta_Lz::Vector{Float64}

    # container_Y # https://jishnub.github.io/SphericalHarmonics.jl/stable/#SphericalHarmonics.cache-Tuple

end

function OrbitalElements_init(nbt::Int64)

    tab_xi = zeros(Float64,2*nbt+1)
    tab_sinhu = zeros(Float64,2*nbt+1)
    tab_tpu = zeros(Float64,2*nbt+1)
    tab_Flmn = zeros(Float64,2*nbt+1)

    tab_alpha_E = zeros(Float64,2*nbt+1)
    tab_alpha_I3 = zeros(Float64,2*nbt+1)
    tab_alpha_Lz = zeros(Float64,2*nbt+1)

    tab_ylm = zeros(Float64,2*nbt+1)
    tab_v = zeros(Float64,2*nbt+1)
    tab_sinv = zeros(Float64,2*nbt+1)
    tab_tpv = zeros(Float64,2*nbt+1)

    tab_beta_E = zeros(Float64,2*nbt+1)
    tab_beta_I3 = zeros(Float64,2*nbt+1)
    tab_beta_Lz = zeros(Float64,2*nbt+1)

    elements = OrbitalElements(0.0,0.0,0.0,0.0,0.0,0.0,0.0,tab_Flmn,tab_xi,tab_sinhu,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_v,tab_sinv,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)

    return elements

end


function OrbitalElements_update!(elem::OrbitalElements, E::Float64, I3::Float64, Lz::Float64, u0::Float64, u1::Float64, v0::Float64, v1::Float64, nbt::Int64, matrix_freq::Matrix{Float64})

    # Ju = _Ju(E,Lz,I3)
    # Jv = _Jv(E,Lz,I3)

    # u0, u1 = find_bounds_u(E,Lz,I3)
    # v0, v1 = find_bounds_v(E,Lz,I3)

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
        # ylm = SphericalHarmonics.associatedLegendre(v,l,m)/sqrt(2.0)
        # ylm = plm/sqrt(2.0)
        # println((ylm1,ylm))

        # println("b2")
        # ylm = elem.container_Y[(l,m)]

        # println("b3")

        # ylm = Ylm(l,m,v,0.0)
        tpv =  _tpv(vEff,E,Lz,I3,v0,v1)

        # println(ylm)
        # elem.tab_ylm[iv+1] = real(ylm)
        # elem.tab_vEff[iv+1] = vEff
        elem.tab_v[iv+1] = v
        elem.tab_sinv[iv+1] = sinv
        # elem.tab_sinvEff[iv+1] = sinvEff

        elem.tab_tpv[iv+1] = tpv

        betaE = Delta^2 *sinv^2/tpv
        betaI3 = Delta^2 /tpv
        betaLz = - Lz/(tpv*sinv^2)

        elem.tab_beta_E[iv+1] = matrix_freq[1,1]*betaE +  matrix_freq[2,1]*betaI3  +  matrix_freq[3,1]*betaLz
        elem.tab_beta_I3[iv+1] = matrix_freq[1,2]*betaE +  matrix_freq[2,2]*betaI3  +  matrix_freq[3,2]*betaLz
        elem.tab_beta_Lz[iv+1] = matrix_freq[1,3]*betaE +  matrix_freq[2,3]*betaI3  +  matrix_freq[3,3]*betaLz

 
    end

    # println("a1")

    for iu=0:2*nbt # reverse sampling
        # println("c0")
        uEff = pi/2 - pi/(2*nbt) * iu
        sinuEff = sin(uEff)
        u = tu*sinuEff + su
        sinhu = sinh(u)
        xi = cosh(u)
        # Flmn = itp_Flmn_xi_exact(l,m,n,xi)
        
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)

        # elem.tab_uEff[iu+1] = uEff
        # elem.tab_u[iu+1] = u 
        elem.tab_xi[iu+1] = xi
        elem.tab_sinhu[iu+1] = sinhu
        # elem.tab_sinuEff[iu+1] = sinuEff

        elem.tab_tpu[iu+1] = tpu

       
        
        # elem.tab_Flmn[iu+1] = Flmn

        # elem.tab_alpha_E[iu+1] = Delta^2 *sinhu^2/tpu
        # elem.tab_alpha_I3[iu+1] = -Delta^2 /tpu
        # elem.tab_alpha_Lz[iu+1] = -Lz/(tpu*sinhu^2)

        alphaE = Delta^2 *sinhu^2/tpu
        alphaI3 = -Delta^2 /tpu
        alphaLz = -Lz/(tpu*sinhu^2)

        elem.tab_alpha_E[iu+1] = matrix_freq[1,1]*alphaE +  matrix_freq[2,1]*alphaI3  +  matrix_freq[3,1]*alphaLz
        elem.tab_alpha_I3[iu+1] = matrix_freq[1,2]*alphaE +  matrix_freq[2,2]*alphaI3  +  matrix_freq[3,2]*alphaLz
        elem.tab_alpha_Lz[iu+1] = matrix_freq[1,3]*alphaE +  matrix_freq[2,3]*alphaI3  +  matrix_freq[3,3]*alphaLz



 
    end

    # println("a2")
    # matrix_freq, matrix_grad = frequency_matrix(E,Lz,I3)

    # elements = OrbitalElements(E,I3,Lz,Ju,Jv,u0,u1,v0,v1,tab_Flmn,tab_uEff,tab_u,tab_xi,tab_sinhu,tab_sinuEff,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_vEff,tab_v,tab_sinv,tab_sinvEff,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)
    # elements = OrbitalElements(E,I3,Lz,u0,u1,v0,v1,tab_Flmn,tab_uEff,tab_u,tab_xi,tab_sinhu,tab_sinuEff,tab_tpu,tab_alpha_E,tab_alpha_I3,tab_alpha_Lz,tab_ylm,tab_vEff,tab_v,tab_sinv,tab_sinvEff,tab_tpv,tab_beta_E,tab_beta_I3,tab_beta_Lz)

    # return elements

end




function OrbitalElements_update_ylm!(elem::OrbitalElements, l::Int64, m::Int64, nbt::Int64)

  
    
    for iv=0:2*nbt # reverse sampling

        v = elem.tab_v[iv+1]
     
        ylm = SphericalHarmonics.associatedLegendre(v,l,m)/sqrt(2.0)
      
        elem.tab_ylm[iv+1] = real(ylm)
 
    end


end


function OrbitalElements_update_Flmn!(elem::OrbitalElements, l::Int64, m::Int64, n::Int64, nbt::Int64)

  
    
    
    for iu=0:2*nbt 
        
        xi = elem.tab_xi[iu+1]
        Flmn = itp_Flmn_xi_exact(l,m,n,xi)
     
        elem.tab_Flmn[iu+1] = Flmn


 
    end


end
