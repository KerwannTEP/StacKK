
println("Nb threads = ",Threads.nthreads())



include("../source/Main.jl")

using HDF5




const nbJ = nbJ_default
const nbt = nbt_default

const Re_omega_max = 1.0#1.0 #1.0
const Re_omega_min = 0.0 # Frequency space symmetrical w.r.t. imaginary axis for non-rotation clusters

const Im_omega_max = 2.0#1.0#2.0
const Im_omega_min = 0.001

const nbRe = 10 #10
const nbIm = 10 #10

const nbgrid = nbRe*nbIm

const tab_absdet_Epq = zeros(Float64, nbgrid)
const tab_argdet_Epq = zeros(Float64, nbgrid)



function tab_det_omega()

    tab_Mpq, tab_omega = ResponseMatrix_m_sampling_alt_half(mmax,Re_omega_min,Re_omega_max,Im_omega_min,Im_omega_max,nbRe,nbIm,nbJ,nbt)

    n0 = 0
    if (mmax != 0)
        n0 = 1
    end

    nbn = nmax - n0 + 1
    nbl = lmax - abs(mmax) + 1

    nbp = nbl * nbn



    for igrid=1:nbgrid


        detEpq = det_Dielectric(tab_Mpq[:,:,igrid], nbp)
        absdet = abs(detEpq)
        argdet = angle(detEpq)

        tab_absdet_Epq[igrid] = absdet
        tab_argdet_Epq[igrid] = argdet
    end


    return tab_absdet_Epq, tab_argdet_Epq, tab_omega


end


# namefile = "../../data/sampling_modes_staeckel_a_"*string(a)*"_m_"*string(mmax)*"_lmax_"*string(lmax)*"_nmax_"*string(nmax)*".hf5"
namefile = "../../data/sampling_modes_staeckel_test_parity_a_"*string(a)*"_m_"*string(mmax)*"_lmax_"*string(lmax)*"_nmax_"*string(nmax)*".hf5"

function writefile!()

    @time tab_absdet_Epq, tab_argdet_Epq, tab_omega = tab_det_omega()

    file = h5open(namefile,"w") # Opening the file

    write(file, "a", a)
    write(file, "nbJ", nbJ)
    write(file, "nbt", nbt)
    write(file, "nbRe", nbRe)
    write(file, "nbIm", nbIm)
    write(file, "re_omega_max", Re_omega_max)
    write(file, "re_omega_min", Re_omega_min)
    write(file, "im_omega_max", Im_omega_max)
    write(file, "im_omega_min", Im_omega_min)

    write(file, "tab_omega_re_im", tab_omega)
    write(file, "tab_absdet", tab_absdet_Epq)
    write(file, "tab_argdet", tab_argdet_Epq)

    write(file, "lmax", lmax)
    write(file, "m", mmax)
    write(file, "nmax", nmax)

    write(file, "kmax", kmax)


    close(file)
end

@time writefile!()
