# include("Args.jl")
# include("Cluster.jl")

using LinearAlgebra

###############################################################
# In-plane motion u 
###############################################################

function Ueff(u::Float64, Lz::Float64, I3::Float64)
    lambda = lambda_from_u(u)

    return Lz^2/(2.0*Delta^2*sinh(u)^4) + (Gtau(lambda)/tanh(u)^2+I3/sinh(u)^2)

end

#  Rewrite this : Simplify MMA ?
function dUeffdu(u::Float64, Lz::Float64, I3::Float64)

    return -4.0*Lz^2*cosh(u)/(2.0*Delta^2*sinh(u)^5) + (sinh(u)*_dUdu(u) - (_U(u)+I3)*2.0*cosh(u))/sinh(u)^3


end

# Minimum energy E_s for an orbit with angular momentum Lz and third integrand I3 
# This corresponds to a shell orbit with constant u=u_s 
# U_s = min_u Ueff(u,Lz,I3)
# Function returns the double (E_s, u_s)
function E_shell(Lz::Float64, I3::Float64)

    b0 = 1.0
    while (dUeffdu(b0,Lz,I3) >= 0.0)
        b0 /= 2.0
    end

    b1 = 1.0
    while (dUeffdu(b1,Lz,I3) <= 0.0)
        b1 *= 2.0
    end
    
    us = bisection(u->dUeffdu(u,Lz,I3),b0,b1)
    Es = Ueff(us,Lz,I3)

    return Es, us 
end



function puSq(u::Float64, E::Float64, Lz::Float64, I3::Float64)

    sinhuSq = sinh(u)^2
    return 2.0*Delta^2*(E*sinhuSq - I3 - _U(u)) - Lz^2/sinhuSq


end

function _tpu(uEff::Float64, E::Float64, Lz::Float64, I3::Float64, u0::Float64, u1::Float64)

    # u0, u1 = find_bounds_u(E,Lz,I3)
    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)
    u = tu*sin(uEff) + su

    # println((puSq(u0,E,Lz,I3),puSq(u1,E,Lz,I3)))

    if (-0.5*PI + cutoffEff < uEff < 0.5*PI - cutoffEff)
        return sqrt(abs(puSq(u,E,Lz,I3)))/(tu*cos(uEff))
    elseif (uEff <= -0.5*PI + cutoffEff)
        uEff_cutoff = -0.5*PI + cutoffEff
        u_cutoff = tu*sin(uEff_cutoff) + su
        tpu_cutoff = sqrt(abs(puSq(u_cutoff,E,Lz,I3)))/(tu*cos(uEff_cutoff))
        dpuSqdu = grad_u_puSq(u0,E,Lz)
        tpu_lim = sqrt(abs(dpuSqdu))/sqrt(abs(u1-u0))
        
        return tpu_lim + (tpu_cutoff-tpu_lim)/cutoffEff * (uEff+0.5*PI)
    else
        uEff_cutoff = 0.5*PI - cutoffEff
        u_cutoff = tu*sin(uEff_cutoff) + su
        tpu_cutoff = sqrt(abs(puSq(u_cutoff,E,Lz,I3)))/(tu*cos(uEff_cutoff))
        dpuSqdu = grad_u_puSq(u1,E,Lz)
        tpu_lim = sqrt(abs(dpuSqdu))/sqrt(abs(u1-u0))
        
        return tpu_cutoff + (tpu_lim-tpu_cutoff)/cutoffEff * (uEff-(0.5*PI - cutoffEff))
    end

end

function grad_u_puSq(u::Float64, E::Float64, Lz::Float64)

    # # don't do the product
    # # simply the cosh as much as possible
    # # to avoid ratio of huge numbers as much as possible
    # # otherwise, we have NaN = Inf*0

    sinhu = sinh(u)
    coshu = cosh(u)

    term0 = 4.0*Delta^2*E*sinhu*coshu
    term1 = 4.0*Delta^2*sinhu*G*M/(c/coshu+sqrt(c^2/coshu^2+Delta^2))
    term2 = Delta^2*G*M*(2.0*Delta^2*sinhu/sqrt(c^2/coshu^2+Delta^2))/(c/coshu+sqrt(c^2/coshu^2+Delta^2))^2
    term3 = 2.0*Lz^2/(sinhu^2*tanh(u))

    return term0 + term1 - term2 + term3

end

# rewrite this for Ueff
# This removes the problematic sinh terms going to infinity
# boundaries are u-solutions to Ueff(u,Lz,I3) = E
# maximum is given by E_shell(Lz,I3)
# special case for Lz=0 ?
function find_bounds_u(E::Float64, Lz::Float64, I3::Float64)

    Es, us = E_shell(Lz,I3)

    # println("sssss:",(E,Es,us,E<Es))
# 
    # println("us=",us)
    # println("Lz/sinh^2=",Lz/sinh(us)^2)


    if (E <= Es) # set non-allowed orbits to circular value by default
        # println("a0")
        # return us, us
        return -1, -1
    else
        # println("b0")
        b1 = us 
        while (Ueff(b1,Lz,I3) < E)
            b1 /= 2.0
        end

        # if (abs(Ueff(us,Lz,I3)-E) <= err_shell) # very close to circular: use circular value
        #     return us, us
        # else

            # println((Ueff(b1,Lz,I3)-E,Ueff(us,Lz,I3)-E))
            u0 = bisection(u->(Ueff(u,Lz,I3)-E),b1,us)

            b2 = us 
            while (Ueff(b2,Lz,I3) < E)
                b2 *= 2.0
            end

            # println((Ueff(us,Lz,I3)-E,Ueff(b2,Lz,I3)-E))

            u1 = bisection(u->(Ueff(u,Lz,I3)-E),us,b2)

            return u0, u1
        # end
    end

end




function _Ju(E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default)

    u0, u1 = find_bounds_u(E,Lz,I3)
    # println("(u0,u1)=",(u0,u1))
    # if (abs(u0-u1) >= err_u) # pv^2(pi/2) >= 0
    if (u0 != -1)# pv^2(pi/2) >= 0
        
        sum = 0.0

        du = (u1-u0)/nbu
        for iu=1:nbu 
            u = u0 + du*(iu-0.5)
            pu = sqrt(abs(puSq(u,E,Lz,I3)))
            sum += pu 
        end

        sum *= du 

        return sum*INV_PI

    else # no orbit
        return 0.0
    end
end


###############################################################
# Out-of-plane motion v 
###############################################################


function pvSq(v::Float64, E::Float64, Lz::Float64, I3::Float64)

    sinvSq = sin(v)^2
    
    return 2.0*Delta^2*(E*sinvSq + I3 + _V(v)) - Lz^2/sinvSq

end


function grad_v_pvSq(v::Float64, E::Float64, Lz::Float64)

    # # don't do the product
    # # simply the cosh as much as possible
    # # to avoid ratio of huge numbers as much as possible
    # # otherwise, we have NaN = Inf*0

    sinv, cosv = sincos(v)

    term0 = 4.0*Delta^2*E*sinv*cosv
    term1 = 4.0*Delta^2*sinv*cosv*G*M/(c+sqrt(c^2+Delta^2*cosv^2))
    term2 = Delta^2*cosv^2*G*M*(-2.0*Delta^2*sinv*cosv/sqrt(c^2+Delta^2*cosv^2))/(c+sqrt(c^2+Delta^2*cosv^2))^2
    term3 = 2.0*cosv*Lz^2/(sinv^3)

    return term0 + term1 + term2 + term3

end


function _tpv(vEff::Float64, E::Float64, Lz::Float64, I3::Float64, v0::Float64, v1::Float64)

    # v0, v1 = find_bounds_v(E,Lz,I3)
    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)
    v = tv*sin(vEff) + sv


    if (-0.5*PI + cutoffEff < vEff < 0.5*PI - cutoffEff)
        return sqrt(abs(pvSq(v,E,Lz,I3)))/(tv*cos(vEff))
    elseif (vEff <= -0.5*PI + cutoffEff)
        vEff_cutoff = -0.5*PI + cutoffEff
        v_cutoff = tv*sin(vEff_cutoff) + sv
        tpv_cutoff = sqrt(abs(pvSq(v_cutoff,E,Lz,I3)))/(tv*cos(vEff_cutoff))
        dpvSqdv = grad_v_pvSq(v0,E,Lz)
        tpv_lim = sqrt(abs(dpvSqdv))/sqrt(abs(v1-v0))
        
        return tpv_lim + (tpv_cutoff-tpv_lim)/cutoffEff * (vEff+0.5*PI)
    else
        vEff_cutoff = 0.5*PI - cutoffEff
        v_cutoff = tv*sin(vEff_cutoff) + sv
        tpv_cutoff = sqrt(abs(pvSq(v_cutoff,E,Lz,I3)))/(tv*cos(vEff_cutoff))
        dpvSqdv = grad_v_pvSq(v1,E,Lz)
        tpv_lim = sqrt(abs(dpvSqdv))/sqrt(abs(v1-v0))
        
        return tpv_cutoff + (tpv_lim-tpv_cutoff)/cutoffEff * (vEff-(0.5*PI - cutoffEff))
    end

end

# special case for Lz=0 ?
function find_bounds_v(E::Float64, Lz::Float64, I3::Float64)

    vm = 0.5*PI
    if (2.0*Delta^2*(E+I3) <= Lz^2) 
        return vm, vm 
    else

        v0 = bisection( v->pvSq(v,E,Lz,I3),0.0,0.5*PI)
        v1 = PI-v0 

        return v0, v1
    end

end

function _Jv(E::Float64, Lz::Float64, I3::Float64, nbv::Int64=nbv_default)

    if ((2.0*Delta^2*(E+I3) >= Lz^2)) # pv^2(pi/2) >= 0
        v0, v1 = find_bounds_v(E,Lz,I3)
        sum = 0.0

        dv = (v1-v0)/nbv
        for iv=1:nbv 
            v = v0 + dv*(iv-0.5)
            pv = sqrt(abs(pvSq(v,E,Lz,I3)))
            sum += pv 
        end

        sum *= dv 

        return sum*INV_PI

    else # no orbit
        return 0.0
    end
end


###############################################################
# Frequencies and gradients
###############################################################

# Semi-log here?
# RK4
function dJudEI3(E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default)
    # dpudE = Delta^2 * sinh(u)^2/pu

    u0, u1 = find_bounds_u(E,Lz,I3)
    


    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    djude = 0.0
    djudi3 = 0.0

 

    for i=1:nbu 
        theta = -0.5*PI + PI/nbu*(i-0.5)
        u = tu*sin(theta) + su
        tpu = _tpu(theta,E,Lz,I3,u0,u1)

        djude += Delta^2 * sinh(u)^2/tpu
        djudi3 += Delta^2/tpu
  
  


    end 

    djude *= 1.0/nbu  # divide by pi
    djudi3 *= -1.0/nbu # divide by pi


    return djude, djudi3
end



# Semi-log
# Do this only for dJudLz
function dJudLz(E::Float64, Lz::Float64, I3::Float64, eps::Float64=epsLz, nbu::Int64=nbu_default)

    u0, u1 = find_bounds_u(E,Lz,I3)

    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    djudlz = 0.0

    # t = log(theta + pi/2)
    # Check that -2.0*Delta^2*(_E0+I3) > 0
    # If not, use tmin_eff = log(pi/(nbu))
    tmin = 0.0
    if (-2.0*Delta^2*(_E0+I3) > 0)
        tmin = log(eps) - log(u1) + 3.0*log(abs(Lz)) - 1.5*log(abs(-2.0*Delta^2*(_E0+I3)))
    else
        tmin = log(PI/nbu)
    end
    tmin_eff = min(tmin,log(PI/nbu)) # Useful for large tmin (I3 large of |Lz| large)
    tmax = log(PI)

    t = tmin_eff
    dt = (tmax-tmin_eff)/nbu

    # Use log sampling on [theta_min, pi/2]
    # RK4: osef ?
    for it=1:nbu 

 
        t = tmin_eff + dt*(it-0.5)
        theta = -0.5*PI + exp(t)
        u = tu*sin(theta) + su
        tpu = _tpu(theta,E,Lz,I3,u0,u1)

        djudlz += -dt*exp(t)*Lz/(tpu*sinh(u)^2)

        
    end

    # Add contribution from [-pi/2 theta_min]
    # mid-point

    # theta_min = exp(tmin)
    dtheta = exp(tmin)

    # Then at theta=-pi/2, use Taylor expansion limit
    theta = -0.5*PI+0.5*dtheta 
    u = tu*sin(theta) + su
    tpu = _tpu(theta,E,Lz,I3,u0,u1)
    djudlz += -dtheta*Lz/(tpu*sinh(u)^2)
 


    djudlz *= 1/pi # divide by pi

    return djudlz
end

function dJu(E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default, eps::Float64=epsLz)

    djude, djudi3 = dJudEI3(E,Lz,I3,nbu)
    djudlz = dJudLz(E,Lz,I3,eps,nbu)

    return djude, djudi3, djudlz
end


# Semi-log here ?
# Not needed I guess?
function dJv(E::Float64, Lz::Float64, I3::Float64, nbv::Int64=nbv_default)


    v0, v1 = find_bounds_v(E,Lz,I3)

    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

    djvde = 0.0
    djvdi3 = 0.0
    djvdlz = 0.0

    for i=1:nbv 
        theta = -0.5*PI + PI/nbv*(i-0.5)
        v = tv*sin(theta) + sv
        tpv = _tpv(theta,E,Lz,I3,v0,v1)

        djvde += Delta^2 * sin(v)^2/tpv
        djvdi3 += Delta^2/tpv
        djvdlz += Lz/(tpv*sin(v)^2)
    end 

    djvde *= 1.0/nbv # divide by pi
    djvdi3 *= 1.0/nbv # divide by pi
    djvdlz *= -1.0/nbv # divide by pi

    return djvde, djvdi3, djvdlz
end


# Semi-log here ?
# Not needed I guess?
function dJvdEI3(E::Float64, Lz::Float64, I3::Float64, nbv::Int64=nbv_default)


    v0, v1 = find_bounds_v(E,Lz,I3)

    sv = 0.5*(v0+v1)
    tv = 0.5*(v1-v0)

    djvde = 0.0
    djvdi3 = 0.0

    for i=1:nbv 
        theta = -0.5*PI + PI/nbv*(i-0.5)
        v = tv*sin(theta) + sv
        tpv = _tpv(theta,E,Lz,I3,v0,v1)

        djvde += Delta^2 * sin(v)^2/tpv
        djvdi3 += Delta^2/tpv
    end 

    djvde *= 1.0/nbv # divide by pi
    djvdi3 *= 1.0/nbv # divide by pi

    return djvde, djvdi3
end

# MANUAL INVERSE 3 By 3 MATRIX 
# IT IS UNSTABLE OTHERWISE ?
# https://ardoris.wordpress.com/2008/07/18/general-formula-for-the-inverse-of-a-3x3-matrix/
function inverse_3_by_3(matrix::Matrix{Float64})

    # @assert (length(matrix[1,:])==3) "Matrix must be 3*3"

    # a, b, c, d, e, f, g, h, i = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    a = matrix[1,1]
    b = matrix[1,2]
    c = matrix[1,3]
    d = matrix[2,1]
    e = matrix[2,2]
    f = matrix[2,3]
    g = matrix[3,1]
    h = matrix[3,2]
    i = matrix[3,3]

    determinant = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)

    inverse = zeros(Float64, 3, 3)

    inverse[1,1] = (e*i - f*h)/determinant
    inverse[1,2] = (c*h - b*i)/determinant
    inverse[1,3] = (b*f - c*e )/determinant
    inverse[2,1] =( f*g - d*i )/determinant
    inverse[2,2] =( a*i - c*g )/determinant
    inverse[2,3] = (c*d - a*f )/determinant
    inverse[3,1] = (d*h - e*g )/determinant
    inverse[3,2] = (b*g - a*h )/determinant
    inverse[3,3] = (a*e - b*d)/determinant

    # inverse /= determinant

    return  inverse 
end


# There might a better way to do this 
# optimize this
function frequency_matrix(E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default, nbv::Int64=nbv_default, eps::Float64=epsLz)

 
        dJudE, dJudI3, dJudLz = dJu(E,Lz,I3,nbu,eps)
        dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3,nbv)



        # matrix_grad = zeros(Float64, 3, 3)

        matrix_grad = zeros(Float64,3,3)
        
        
        matrix_grad[1,1] = dJudE 
        matrix_grad[1,2] = dJudI3 
        matrix_grad[1,3] = dJudLz
        matrix_grad[2,1] = dJvdE 
        matrix_grad[2,2] = dJvdI3 
        matrix_grad[2,3] = dJvdLz
        matrix_grad[3,3] = 1.0



        return inverse_3_by_3(matrix_grad), matrix_grad

end


function fill_grad_frequency_matrix!(matrix_grad::Matrix{Float64}, matrix_freq::Matrix{Float64},  E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default, nbv::Int64=nbv_default, eps::Float64=epsLz)

 
    dJudE, dJudI3, dJudLz = dJu(E,Lz,I3,nbu,eps)
    dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3,nbv)

    # matrix_grad = zeros(Float64, 3, 3)    
    
   

    # a = dJudE 
    # b = dJudI3 
    # c = dJudLz
    # d = dJvdE 
    # e = dJvdI3 
    # f = dJvdLz
    # g = 0.0
    # h = 0.0
    # i = 1.0

    matrix_grad[1,1] = dJudE
    matrix_grad[1,2] = dJudI3 
    matrix_grad[1,3] = dJudLz
    matrix_grad[2,1] = dJvdE 
    matrix_grad[2,2] = dJvdI3 
    matrix_grad[2,3] = dJvdLz
    matrix_grad[3,1] = 0.0
    matrix_grad[3,2] = 0.0
    matrix_grad[3,3] = 1.0


    determinant = dJudE *dJvdI3  - dJudI3 *dJvdE 

    # inverse = zeros(Float64, 3, 3)

    matrix_freq[1,1] = dJvdI3/determinant
    matrix_freq[1,2] = -dJudI3 /determinant
    matrix_freq[1,3] = (dJudI3 *dJvdLz - dJudLz*dJvdI3  )/determinant
    matrix_freq[2,1] = -dJvdE/determinant
    matrix_freq[2,2] =  dJudE /determinant
    matrix_freq[2,3] = (dJudLz*dJvdE  - dJudE * dJvdLz )/determinant
    matrix_freq[3,1] = 0.0
    matrix_freq[3,2] = 0.0
    matrix_freq[3,3] = 1.0


    

end

function frequency_matrix_test(E::Float64, Lz::Float64, I3::Float64, nbu::Int64=nbu_default, nbv::Int64=100, eps::Float64=0.01)

     dJudE, dJudI3, dJudLz = dJu(E,Lz,I3,nbu,eps)
     dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3,nbv)

    detJ = dJudE*dJvdI3-dJudI3*dJvdE

    Omega_u = dJvdI3/detJ 
    Omega_v = -dJudI3/detJ 
    Omega_z = (dJudI3*dJvdLz-dJudLz*dJvdI3)/detJ 

    dI3dJu = -dJvdE/detJ 
    dI3dJv = dJudE/detJ 
    dI3dLz = (dJudLz*dJvdE-dJudE*dJvdLz)/detJ

    println((Omega_u, Omega_v, Omega_z))


    return Omega_u, Omega_v, Omega_z, dI3dJu, dI3dJv, dI3dLz

end



# Ubar(u,Lz) = U(u) + Lz^2/(2 Delta^2 sinh^2 u
function Ubar(u::Float64, Lz::Float64)

    return cosh(u)^2*(-G*M)/(c + sqrt(c^2 + (a^2 - c^2) *cosh(u)^2)) + 
    Lz^2/(2.0 * (a^2 - c^2) *sinh(u)^2)

end

function d2Ubardu2(u::Float64, Lz::Float64)

    num1 = G *M* (4.0 *(a^2 + c^2) *cosh(2 *u) + (a - c)* (a + c) *(3.0 + cosh(4.0 *u)))
    den1 = 2.0* sqrt(2.0)* (a^2 + c^2 + (a - c)* (a + c) *cosh(2* u))^(1.5)

    num2 = Lz^2* (2.0 + cosh(2.0 *u))* 1.0/sinh(u)^4
    den2 = a^2 - c^2


    return -num1/den1 + num2/den2


end

###############################################################
# Utilitary functions 
###############################################################


function bisection(fun::Function, xl::Float64, xu::Float64, tolx::Float64=1.0*10^(-10), tolf::Float64=1.0*10^(-10), iterMAX::Int64=50)
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####

    # if (abs(fl) <= tolf) # We have already found a solution on the left bracket
    #     return xl # Returning the left bracket
    # end
    # #####
    # if (abs(fu) <= tolf) # We have already found a solution on the right bracket
    #     return xu # Returning the right bracket
    # end

    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET : (xl,xu,fl,fu) = "*string((xl,xu,fl,fu))
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if ((abs(xu-xl) <= tolx) || (iter > iterMAX)) # The considered bracket is smaller than the tolerance, or we have made too many iterations
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the counter of iterations
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
end


# Backward integration
# USE RK4
# Tqke inspiration from https://github.com/JuliaStellarDynamics/OrbitalElements.jl/blob/main/src/Utils/Integrators.jl
function tab_alphak(k1::Int64, k2::Int64, k3::Int64, Ju::Float64, Jv::Float64, Lz::Float64, nbu::Int64=100)

    E, Lz, I3 = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    u0, u1 = find_bounds_u(E,Lz,I3)
    su = 0.5*(u0+u1)
    tu = 0.5*(u1-u0)

    # println((u0,u1))

    trans_freq = transpose(frequency_matrix(E,Lz,I3)[1])
    trans_bk = [k1 k2 k3]

    # display(trans_freq)

    sum = [0.0, 0.0, 0.0]

    tab_alphak = zeros(Float64, 2, nbu) # (u, alphak[u] )

    uEff = pi/2
    du = pi/nbu

    # Test 

    tabint = zeros(Float64, 2, nbu)
    


    for i=1:nbu 
        # uEff = pi/2 - pi/nbu*(i-0.5)

        # Step 1
        u = tu*sin(uEff) + su
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)
        
        sum[1] += du/6.0 * Delta^2 *sinh(u)^2/tpu
        sum[2] += -du/6.0 * Delta^2 /tpu
        sum[3] += -du/6.0 * Lz/(tpu*sinh(u)^2)

        # Step 2 and 3
        uEff -= du/2.0
        u = tu*sin(uEff) + su
        tpu = _tpu(uEff,E,Lz,I3,u0,u1)

        tabint[1,i],  tabint[2,i] = uEff+pi/2, Lz/(tpu*sinh(u)^2)

        sum[1] += 4.0*du/6.0 * Delta^2 *sinh(u)^2/tpu
        sum[2] += -4.0*du/6.0 * Delta^2 /tpu
        sum[3] += -4.0*du/6.0 * Lz/(tpu*sinh(u)^2)

         # Step 4
         uEff -= du/2.0
         u = tu*sin(uEff) + su
         tpu = _tpu(uEff,E,Lz,I3,u0,u1)
 
         sum[1] += du/6.0 * Delta^2 *sinh(u)^2/tpu
         sum[2] += -du/6.0 * Delta^2 /tpu
         sum[3] += -du/6.0 * Lz/(tpu*sinh(u)^2)


        # println((sum[1],sum[2],sum[3]))


 
        # display(trans_freq*sum)
        # display(bk)

        akr = trans_bk*trans_freq*sum
        
        tab_alphak[1,i], tab_alphak[2,i] = uEff, k1*pi-akr[1]
    end

    # pt = plot(tabint[1,:],  tabint[2,:], xaxis=:log10)
    # savefig(pt,"test_integrand.png")

    return tab_alphak

end
