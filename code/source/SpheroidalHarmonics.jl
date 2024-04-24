using SphericalHarmonics


# https://juliapackages.com/p/sphericalharmonics
function Ylm(ell::Int64, m::Int64, v::Float64, phi::Float64)

    # eta = cos(v)
    # Y = computeYlm(v, phi, lmax = ell)
    Y = computeYlm(v, phi, ell, m)
    return Y[(ell,m)]
end


