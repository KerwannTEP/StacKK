# include("Args.jl")
# include("Constants.jl")

function _tE(E::Float64)

    return -(a+c)*E/(G*M)
end 

function _tLz(Lz::Float64)

    return Lz/sqrt((a+c)*G*M)
end 

