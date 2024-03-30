using HDF5

# include("SpheroidalHarmonics.jl")

# Read the radial orthogonalized basis elements from the .hdf5 files
# Interpolate linearly the functions
# Radial basis function are orthogonalized for the scalar product
#
# <F^lmk, F^lmn> = - int \rd \xi F^lmk D^lmn
#
# The global potential basis functions are given by
#
# psi^lmn = sqrt(4 pi G)/Delta F^lmn(\xi) Y_l^m(v,\phi)
#
# The global scalar product is given by
#
# <psi^p, psi^q> = int \rd \br \psi^p(\br)^{*} \rho^q(\br)
#
# Using the F^lmn orthogonalized as in described, the total scalar product yields
#
# <psi^p, psi^q> = 1/Delta delta_p^q
#
#
#
# This is slightly difference from the usual radial scalar product: <F^lmk, F^lmn> = - 1/Delta int \rd \xi F^lmk D^lmn
#
# This one yields the total scalar product
#
# <psi^p, psi^q> = delta_p^q

# Get the data file
function get_tab_Flmn_exact(l::Int64, m::Int64, n::Int64)


    namefile = "../../basis_functions_h_1/l_"*string(l)*"/m_"*string(abs(m))*"/F_l_"*string(l)*"_m_"*string(m)*"_n_"*string(n)*".hdf5"

    file = h5open(namefile)
    data = read(file, "InterpolationTable")

    # x = (xi-1)/(xi+1)
    data_x = data[1,:]
    data_F = data[2,:]

    return data_x, data_F
end


function get_all_tab_Flmn_exact()

    tab_x = Vector{Float64}[]
    tab_F = Vector{Float64}[]

    tab_lmn = Any[]
    tab_p = zeros(Int64, lmax+1, mmax+1, nmax+1)

    nb_lmn = 0
    index_p = 1

    println("Loading radial functions tables")

    for l=0:lmax
        for m=0:min(l,mmax)
            n0 = 0
            if (m > 0)
                n0 = 1
            end
            for n=n0:nmax

                push!(tab_lmn,[l,m,n])

                tab_p[l+1,m+1,n+1] = index_p

                index_p += 1

                nb_lmn += 1

                
            end

        end
    end

    for p=1:nb_lmn

        l, m , n = tab_lmn[p]

        println("(l,m,n) = ",(l,m,n))


        namefile = "../../basis_functions_h_1/l_"*string(l)*"/m_"*string(abs(m))*"/F_l_"*string(l)*"_m_"*string(m)*"_n_"*string(n)*".hdf5"
        file = h5open(namefile)
        data = read(file, "InterpolationTable")

        push!(tab_x,data[1,:])
        push!(tab_F,data[2,:])

    end


    return tab_x, tab_F, tab_lmn, tab_p

end

const tab_x, tab_F, tab_lmn, tab_p = get_all_tab_Flmn_exact()

# Linear interpolation of Flmn(x)
# x = (xi-1)/(xi+1) in [0, 1]
# xi in [1, infty]
# function itp_Flmn_x_exact(x::Float64, data_x, data_F)

#     nbxi = length(data_x)

#     # list starts at 0 and ends at 1: nbxi points, hence nbxi-1 intervals
#     # dx = 1.0/(nbxi-1)
#     ix = ceil(Int64,x*(nbxi-1)) # Left index

#     if (x < 1.0)
#         xleft = data_x[ix]
#         xright = data_x[ix+1]
#         Fleft = data_F[ix]
#         Fright = data_F[ix+1]

#         slope = (Fright-Fleft)/(xright-xleft)

#         return sqrt(Delta)*(Fleft + (x-xleft)*slope )
#     else
#         return 0.0
#     end
# end

# function finite_diff(p::Int64, ix::Int64, x::Float64)

#     xleft = tab_x[p][ix]
#     xright = tab_x[p][ix+1]
#     Fleft = tab_F[p][ix]
#     Fright = tab_F[p][ix+1]

#     slope = (Fright-Fleft)/(xright-xleft)

#     Fx =  sqrtDelta*(Fleft + (x-xleft)*slope )
    
#     return Fx

# end



function itp_Flmn_x_exact(l::Int64, m::Int64, n::Int64, x::Float64)


    p = tab_p[l+1,abs(m)+1,n+1]
    nbxi = n*100 # length(tab_x[p])
    if (n==0)
        nbxi = 100
    end

    # tabx = tab_x[p]
    # tabF = tab_F[p]

    # list starts at 0 and ends at 1: nbxi points, hence nbxi-1 intervals
    # dx = 1.0/(nbxi-1)
    ix = ceil(Int64,x*(nbxi-1)) # Left index

    if (x < 1.0)
        xleft = tab_x[p][ix]
        xright = tab_x[p][ix+1]
        Fleft = tab_F[p][ix]
        Fright = tab_F[p][ix+1]

        slope = (Fright-Fleft)/(xright-xleft)

        Fx =  sqrtDelta*(Fleft + (x-xleft)*slope )
    
    return Fx
       
    else
        return 0.0
    end
end

# x = (xi-1)/(xi+1)
# Orthonormalization requires to multiply by sqrt(Delta)
# This is because the file contains radial functions which obeys
# < F^lmk, F^lmn> = 1/Delta delta_k^n
# function itp_Flmn_xi_exact(xi::Float64, data_x, data_F)

#     x = (xi-1)/(xi+1)

#     return itp_Flmn_x_exact(x,data_x,data_F)
# end

function itp_Flmn_xi_exact(l::Int64, m::Int64, n::Int64, xi::Float64)

    x = (xi-1)/(xi+1)

    return itp_Flmn_x_exact(l,m,n,x)
end

# function psi_lmn_exact(xi::Float64, eta::Float64, phi::Float64, l::Int64, m::Int64, n::Int64, data_x, data_F)

#     Flmn = itp_Flmn_xi_exact(xi,data_x,data_F)

#     v = acos(eta)
#     Ylm = Ylm(l,m,v,phi)

#     return sqrt(4.0*pi*G)/Delta*Flmn*Ylm
# end
