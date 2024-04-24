# include("AnglesActions.jl")
# include("Constants.jl")
# include("Args.jl")

# brackets errors for E_shell sometimes
function E_Lz_I3_from_Ju_Lz_Jv(Ju::Float64, Lz::Float64, Jv::Float64, nbu::Int64=nbu_default, err::Float64=errInverse, maxIter::Int64=maxIter_default)

    # Initialize within allowed (E, I3) region
    I3_guess = Lz^2/(2.0*Delta^2) - _E0 # E + I3 = E - E_0 + Lz^2/(2*Delta^2) >= Lz^2/(2*Delta^2)
    E_guess = E_shell(Lz,I3_guess)[1]/2.0

    Ju_guess = _Ju(E_guess,Lz,I3_guess,nbu)
    Jv_guess = _Jv(E_guess,Lz,I3_guess,nbu)

    iter = 0

    # Newton's method
    while (((Ju_guess-Ju)^2 + (Jv_guess-Jv)^2 > err^2) && (iter < maxIter))

        dJudE, dJudI3 = dJudEI3(E_guess,Lz,I3_guess,nbu)
        dJvdE, dJvdI3 = dJvdEI3(E_guess,Lz,I3_guess,nbu)

        determinant = dJudE *dJvdI3  - dJudI3 *dJvdE 

        # Manually compute the inverse to save allocations and avoid huge garbage collection (GC)
        invJacobian11 = dJvdI3/determinant
        invJacobian12 = -dJudI3 /determinant
        invJacobian21 = -dJvdE/determinant
        invJacobian22 =  dJudE /determinant     

        delJu = Ju_guess-Ju
        delJv = Jv_guess-Jv
        
        next_E = E_guess - (invJacobian11 * delJu + invJacobian12 * delJv )
        next_I3 = I3_guess - (invJacobian21 * delJu + invJacobian22 * delJv )

        # Do not exit the allow (E,Lz,I3) space
        while (!((next_E >= E_shell(Lz,next_I3)[1]) && (next_E+next_I3>=Lz^2/(2.0*Delta^2))))
            next_E = 0.5*(next_E+E_guess)
            next_I3 = 0.5*(next_I3+I3_guess)
        end

        # Do not exit the allow (E,Lz,I3) space
        while (next_E >= 0)
            next_E = 0.5*(next_E+E_guess)
            next_I3 = 0.5*(next_I3+I3_guess)
        end

        # Update guess
        E_guess = next_E
        I3_guess = next_I3

        Ju_guess = _Ju(E_guess,Lz,I3_guess,nbu)
        Jv_guess = _Jv(E_guess,Lz,I3_guess,nbu)

        iter += 1

    end

    return E_guess, Lz, I3_guess
end


##########################################
# Frequencies from actions 
##########################################

# Frequencies are in the first line of the inverse of the Jacobian matrix of the mapping (E,Iz,Lz) -> (Ju,Jv,Lz)
# Returns (Omega_u, Omega_v, Omega_phi)
function frequency_vector(Ju::Float64, Jv::Float64, Lz::Float64)

    # println("a0")
    E, Lz, I3, _ = E_Lz_I3_from_Ju_Lz_Jv(Ju,Lz,Jv)

    dJudE, dJudI3, dJudLz = dJu(E,Lz,I3)
    dJvdE, dJvdI3, dJvdLz = dJv(E,Lz,I3)

    Jacobian = [dJudE dJudI3 dJudLz; dJvdE dJvdI3 dJvdLz; 0 0 1]
    invJacobian = inv(Jacobian)

    display(Jacobian)

    return invJacobian[1,1], invJacobian[1,2], invJacobian[1,3]
end


