mutable struct OrbitalElements_alt_half

    E::Float64
    I3::Float64 
    Lz::Float64

    u0::Float64 
    u1::Float64 

    v0::Float64 
    v1::Float64

    tab_Flmn::Vector{Float64}
    tab_tpu::Vector{Float64}
    tab_dpudE::Vector{Float64}
    tab_dpudI3::Vector{Float64}
    tab_dpudLz::Vector{Float64}
    tab_Jacu::Vector{Float64}

    tab_ylm::Vector{Float64}
    tab_tpv::Vector{Float64}
    tab_dpvdE::Vector{Float64}
    tab_dpvdI3::Vector{Float64}
    tab_dpvdLz::Vector{Float64}
    tab_Jacv::Vector{Float64}


end


function OrbitalElements_alt_half_init(nbt::Int64)

    tab_Flmn = zeros(Float64, 2*nbt+1)
    tab_tpu = zeros(Float64, 2*nbt+1)
    tab_dpudE = zeros(Float64, 2*nbt+1)
    tab_dpudI3 = zeros(Float64, 2*nbt+1)
    tab_dpudLz = zeros(Float64, 2*nbt+1)
    tab_Jacu = zeros(Float64, 2*nbt+1)

    tab_ylm = zeros(Float64, 2*nbt+1)
    tab_tpv = zeros(Float64, 2*nbt+1)
    tab_dpvdE = zeros(Float64, 2*nbt+1)
    tab_dpvdI3 = zeros(Float64, 2*nbt+1)
    tab_dpvdLz = zeros(Float64, 2*nbt+1)
    tab_Jacv = zeros(Float64, 2*nbt+1)

    

    elements = OrbitalElements_alt_half(0.0,0.0,0.0,0.0,0.0,0.0,0.0,tab_Flmn,tab_tpu,tab_dpudE,tab_dpudI3,tab_dpudLz,tab_Jacu,tab_ylm,tab_tpv,tab_dpvdE,tab_dpvdI3,tab_dpvdLz,tab_Jacv)

    return elements

end


function OrbitalElements_alt_half_update!(elem::OrbitalElements_alt_half, l::Int64, m::Int64, n::Int64, E::Float64, I3::Float64, Lz::Float64, u0::Float64, u1::Float64, v0::Float64, v1::Float64, grad_matrix::Matrix{Float64}, nbt::Int64)

    elem.E = E 
    elem.Lz = Lz 
    elem.I3 = I3

    elem.u0 = u0
    elem.u1 = u1
    elem.v0 = v0
    elem.v1 = v1

    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

    Jw = abs(grad_matrix[1,1]*grad_matrix[2,2] - grad_matrix[1,2]*grad_matrix[2,1])


    for iu=0:2*nbt

        uEff = pi/2 - pi/(2.0*nbt)*(iu)
        u = tu*sin(uEff) + su
        xi = cosh(u)
        Flmn_u = Flmn(l,m,n,xi)#itp_Flmn_xi_exact(l,m,n,xi)
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)

        dpudE = Delta^2*sinh(u)^2/tpu
        dpudI3 = -Delta^2/tpu
        dpudLz = -Lz/(tpu*sinh(u)^2)
        Jacu = Delta^4/Jw * (sinh(u)^2 ) 

        elem.tab_Flmn[iu+1] = Flmn_u
        elem.tab_tpu[iu+1] = tpu
        elem.tab_dpudE[iu+1] = dpudE
        elem.tab_dpudI3[iu+1] = dpudI3
        elem.tab_dpudLz[iu+1] = dpudLz
        elem.tab_Jacu[iu+1] = Jacu


    end

    for iv=0:2*nbt

        vEff = 0.0 - pi/2.0/(2.0*nbt)*(iv)
        v = tv*sin(vEff) + sv
        ylm_v = SphericalHarmonics.associatedLegendre(v,l,m)/sqrt(2.0)
        tpv = _tpv(vEff,E,Lz,I3,v0,v1)
        
        dpvdE = Delta^2*sin(v)^2/tpv
        dpvdI3 = Delta^2/tpv 
        dpvdLz = -Lz/(tpv*sin(v)^2)
        Jacv = Delta^4/Jw * ( sin(v)^2) 


        elem.tab_ylm[iv+1] = ylm_v
        elem.tab_tpv[iv+1] = tpv
        elem.tab_dpvdE[iv+1] = dpvdE
        elem.tab_dpvdI3[iv+1] = dpvdI3
        elem.tab_dpvdLz[iv+1] = dpvdLz
        elem.tab_Jacv[iv+1] = Jacv


    end

   

end