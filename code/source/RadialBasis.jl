using HDF5 

# include("SpheroidalHarmonics.jl")

# Read the radial orthogonalized basis elements from the .hdf5 files 
# Interpolate linearly the functions 

# Get the data file
function get_tab_Flmn(l::Int64, m::Int64, n::Int64)

    
    namefile = "../../nb/potential/h_1/l_"*string(l)*"/m_"*string(m)*"/F_l_"*string(l)*"_m_"*string(m)*"_n_"*string(n)*".hdf5"

    file = h5open(namefile)
    data = read(file, "InterpolationTable")

    # x = (xi-1)/(xi+1)
    data_x = data[1,:]
    data_F = data[2,:]

    return data_x, data_F 
end



# Linear interpolation of Flmn(x)
# x = (xi-1)/(xi+1) in [0, 1]
# xi in [1, infty]
function itp_Flmn_x(x::Float64, data_x, data_F)

    nbxi = length(data_x)

    # list starts at 0 and ends at 1: nbxi points, hence nbxi-1 intervals
    # dx = 1.0/(nbxi-1) 
    ix = ceil(Int64,x*(nbxi-1)) # Left index 

    if (x < 1.0)
        xleft = data_x[ix]
        xright = data_x[ix+1]
        Fleft = data_F[ix]
        Fright = data_F[ix+1]

        slope = (Fright-Fleft)/(xright-xleft)

        return Fleft + (x-xleft)*slope 
    else
        return 0.0
    end
end

# x = (xi-1)/(xi+1)
function itp_Flmn_xi(xi::Float64, data_x, data_F)

    x = (xi-1)/(xi+1)

    return itp_Flmn_x(x,data_x,data_F)
end

function psi_lmn(xi::Float64, eta::Float64, phi::Float64, l::Int64, m::Int64, n::Int64, data_x, data_F)

    Flmn = itp_Flmn_xi(xi,data_x,data_F)

    v = acos(eta)
    Ylm = Ylm(l,m,v,phi)

    return sqrt(4.0*pi*G)/Delta*Flmn*Ylm
end