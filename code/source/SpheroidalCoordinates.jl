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



function test_coordinates(nbx::Int64, nbeta::Int64, nbphi::Int64)

    # Integrate rho(r)

    # R,z,phi 

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

    println("spheroidal = ",sum)
end



