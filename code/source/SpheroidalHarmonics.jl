using SphericalHarmonics


# https://juliapackages.com/p/sphericalharmonics
function Ylm(ell::Int64, m::Int64, v::Float64, phi::Float64)

    # eta = cos(v)
    # Y = computeYlm(v, phi, lmax = ell)
    Y = computeYlm(v, phi, ell, m)
    return Y[(ell,m)]
end





# Table of Ylm(l,m,v,phi=0) for l=0,1,...,lmax and m=0,1,...,lmax
# Ylm(l,-m,v,0.0) = (-1)^m Ylm(l,m,v,0.0) 
# No need for conjugate since phi=0 yields real spherical harmonics
# Ylm = 0.0 for m>l

# const tab_Ylm = zeros(Float64, lmax+1, lmax+1, nbvInt+1)  

# function tab_Ylm!()

#     for l=0:lmax
#         for m=0:l
#             for iv=1:nbvInt+1
#                 v = pi/(nbvInt) * (iv-1)
#                 ylm = Ylm(l,m,v,0.0)
#                 tab_Ylm[l+1,m+1,iv] = real(ylm)
#             end
#         end
#     end


# end

# @time tab_Ylm!()

# # Interpolate Ylm(l,m,v,0.0)
# function itp_ylm(l::Int64, m::Int64, v::Float64)

#     dv = pi/(nbvInt)
#     iv = ceil(Int64,v/dv) 

#     abs_m = abs(m)
#     sgn_m = sign(m)


#     ylm_inf = tab_Ylm[l+1,abs_m+1,iv]
#     ylm_sup = tab_Ylm[l+1,abs_m+1,iv+1]

#     slope = (ylm_sup-ylm_inf)/(dv)

#     ylm = ylm_inf + (v-dv*iv)*slope 

#     if (sgn_m > 0)
#         return ylm 
#     else
#         return (-1)^m * ylm 
#     end
# end


