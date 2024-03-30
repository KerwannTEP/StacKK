include("../source/Main.jl")

using HDF5


const nbJ = 10
const nbt = 100

const Re_omega_max = 0.5
const Re_omega_min = -0.5

const Im_omega_max = 0.5
const Im_omega_min = 0.01

const nbRe = 2
const nbIm = 2

const nbgrid = nbRe*nbIm

const tab_omega_re_im = zeros(Float64, 2, nbgrid)
const tab_absdet_Epq = zeros(Float64, nbgrid)

function tab_omega_re_im!()

    igrid = 1
    for ire=1:nbRe 
        re_omega = Re_omega_min + (Re_omega_max-Re_omega_min)/(nbRe-1)*(ire-1)
        for iim=1:nbIm
            im_omega = Im_omega_min + (Im_omega_max-Im_omega_min)/(nbIm-1)*(iim-1)

            tab_omega_re_im[1,igrid] = re_omega
            tab_omega_re_im[2,igrid] = im_omega
            
        
        end
    end

end

@time tab_omega_re_im!()

function tab_absdet_Epq!()

    for igrid=1:nbgrid 

        re_omega = tab_omega_re_im[1,igrid] 
        im_omega = tab_omega_re_im[2,igrid]

        omega = re_omega + 1im*im_omega

        tab_Mpq, tab_ln = ResponseMatrix_m(mmax,omega,nbJ,nbt)
        size = length(tab_ln[:,1])
        absdet = abs(det_Dielectric(tab_Mpq, size))

        tab_absdet_Epq[igrid] = absdet 
    end




    
end 

@time tab_absdet_Epq!()

namefile = "../../data/sampling_modes_staeckel_a_"*string(a)*"_m_"*string(mmax)*"_lmax_"*string(lmax)*"_nmax_"*string(nmax)*".hf5"

function writefile!()

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

    write(file, "tab_omega", tab_omega_re_im)
    write(file, "tab_absdet", tab_absdet_Epq)

    write(file, "lmax", lmax)
    write(file, "m", mmax)
    write(file, "nmax", nmax)

    write(file, "kmax", kmax)
    

    close(file)
end