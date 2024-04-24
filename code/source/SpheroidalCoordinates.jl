# include("Args.jl")

function R_z_from_lambda_nu(lambda::Float64, nu::Float64)

    RSq = (lambda+alpha)*(nu+alpha)/(alpha-gamma)
    zSq = (lambda+gamma)*(nu+gamma)/(gamma-alpha)

    return sqrt(RSq), sqrt(zSq)
end 

function lambda_nu_from_u_v(u::Float64, v::Float64)

    lambda = a^2*cosh(u)^2 - c^2*sinh(u)^2
    nu = c^2*sin(v)^2 + a^2*cos(v)^2

    return lambda, nu 
end 

function lambda_from_u(u::Float64)

    # lambda = a^2*cosh(u)^2 - c^2*sinh(u)^2
    lambda = a^2+ Delta^2*sinh(u)^2


    return lambda
end 

function nu_from_v(v::Float64)


    nu = c^2*sin(v)^2 + a^2*cos(v)^2

    return nu 
end 

function R_z_from_u_v(u::Float64, v::Float64)

    R = Delta*sinh(u)*sin(v)
    z = Delta*cosh(u)*cos(v)

    return R, z
end

function xi_eta_from_u_v(u::Float64, v::Float64)

    xi = cosh(u)
    eta = cos(v)

    return xi, eta
end






