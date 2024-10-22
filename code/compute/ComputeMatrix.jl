# Compute the response matrix at a given frequency 
println("Nb threads = ",Threads.nthreads())

include("../source/Main.jl")

const omega = 0.0 + 0.2im 

# WRITE A MATRIX EVALUATION WITHOUT GRID 
# ADAPT GRID FUNCTION FOR 1*1 GRID CASE 
# COMPARE WITH TEX_STACKEL


# tab_Mpq, tab_omega = ResponseMatrix_m_GS_sampling(mmax, real(omega), real(omega), imag(omega), imag(omega), 1, 1, nbJ_default, nbt_default)
tab_Mpq_a, tab_Mpq_b, tab_Mpq_c, tab_omega = ResponseMatrix_m_GS_sampling_rot_split_separate(mmax,  real(omega), real(omega), imag(omega), imag(omega), 1, 1, nbJ_default, nbJ_default, nbJ_default, 
                                                                                            nbt_default, 1.0*10^(-3), 1, nbJ_default^3, 1, nbJ_default^2)
# display(tab_Mpq)

display(tab_Mpq_a)
display(tab_Mpq_b)
display(tab_Mpq_c)
