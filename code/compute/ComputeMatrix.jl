# Compute the response matrix at a given frequency 


println("Nb threads = ",Threads.nthreads())

include("../source/Main.jl")

const omega = 0.0 + 0.2im 
const alphaRot = 1.0

# Get the GS decomposition out of the matrix
# Use only for nmax <= 10
# Beyond that: numerical instability
tab_Mpq_a, tab_Mpq_b, tab_Mpq_c = ResponseMatrix_m_GS_rot_separate(mmax, omega, nbJ_default, nbJ_default, nbJ_default, nbt_default, 1.0*10^(-3))
tab_Mpq = tab_Mpq_a + alphaRot * (tab_Mpq_b + tab_Mpq_c)

display(tab_Mpq)

println("---------------------")

# Use the GS decomposition within the matrix 
# Slower but numerically accurate
tab_Mpq_a, tab_Mpq_b, tab_Mpq_c = ResponseMatrix_m_rot_separate(mmax, omega, nbJ_default, nbJ_default, nbJ_default, nbt_default, 1.0*10^(-3))
tab_Mpq = tab_Mpq_a + alphaRot * (tab_Mpq_b + tab_Mpq_c)

display(tab_Mpq)