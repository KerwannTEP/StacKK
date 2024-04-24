# include("Args.jl")
# include("Constants.jl")
# include("SpheroidalCoordinates.jl")

function Gtau(tau::Float64)

    
    return -G*M/(c+sqrt(tau))
end

function GtauGrad(tau::Float64)

    return G*M*(0.5/sqrt(tau))/(c+sqrt(tau))^2
end

function psi(lambda::Float64, nu::Float64)
    
    return -G*M/(sqrt(lambda)+sqrt(nu))
end 

function rho(lambda::Float64, nu::Float64)

    num = lambda*nu + a^2*(lambda+3.0*sqrt(lambda*nu)+nu)
    den = (lambda*nu)^(3/2)*(sqrt(lambda)+sqrt(nu))^3

    return M*c^2/(4.0*pi) * num/den 
end


function psi_u_v(u::Float64, v::Float64)

    lambda, nu = lambda_nu_from_u_v(u,v)

    return psi(lambda,nu)

end

function _U(u::Float64)

    lambda = lambda_from_u(u)

    return cosh(u)^2*Gtau(lambda)

end

function _dUdu(u::Float64)

    lambda = lambda_from_u(u)
    dlambdadu = 2.0*Delta^2*cosh(u)*sinh(u)

    return 2.0*cosh(u)*sinh(u)*Gtau(lambda) + cosh(u)^2*dlambdadu*GtauGrad(lambda)
end

function _V(v::Float64)

    nu = nu_from_v(v)

    return cos(v)^2 * Gtau(nu)

end
