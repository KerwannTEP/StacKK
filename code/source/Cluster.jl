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

function psi_u_v_test(u::Float64, v::Float64)

    return (_U(u)-_V(v))/(sinh(u)^2+sin(v)^2)
end

# integrate rho(xi,eta)

function integrate_rho(nbx::Int64=100, nbeta::Int64=100)

    sum = 0.0

    # int \rd \br rho(\br) = \Delta int \rd \xi \rd \eta \rd \phi c(\xi,\eta) \rho(\xi,\eta)
    # x = (xi-1)/(xi+1) in [0, 1]
    # xi = (1+x)/(1-x)
    # dxi = 2/(1-x)^2 dx

    # c(xi,eta) = Delta^2 xi^2 - Delta^2 eta^2

    for i=1:nbx 
        x = 1.0/nbx*(i-0.5)
        xi = (1+x)/(1-x)
        lambda = c^2 + Delta^2*xi^2

        sumx = 0.0

        for j=1:nbeta 
            eta = -1.0 + 2.0/nbeta*(j-0.5)
            nu = c^2 + Delta^2*eta^2
            v = acos(eta)

            rhor = rho(lambda,nu)
            cxieta = Delta^2*(xi^2 - eta^2)

            # println((rhor,cxieta,2.0/(1-x)^2))

            sumx += 2.0/nbeta * cxieta * rhor
        end 

        sum += 1.0/nbx * 2.0/(1-x)^2 * sumx 
    end

    sum *= 2.0*pi # integration over \phi yields 2*pi

    sum *= Delta

    return sum 
end
