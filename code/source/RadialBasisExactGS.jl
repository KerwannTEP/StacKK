using HDF5

function Flmn_bare(xi::Float64, l::Int64, m::Int64, n::Int64)

    h = 1.0
    Wrad = 1.0/(xi+1)
    flmn = 1.0
    if (n>0)
        flmn = ((xi-1.0)/(xi+h))^n
    end

    return Wrad*flmn 
end


function C0n(n::Int64)
    h=1.0
    return 2.0*(h+1.0)^2*n^2
end

function C1n(n::Int64)
    h=1.0
    p=1.0
    return (1.0+(h+1.0)*n*(n+1.0) - (2.0*n+1.0)*(2.0*p+1.0))*(h+1.0)
end

function C2n(n::Int64)
    h=1.0
    p=1.0
    return 2.0*p*(p-(h+1.0)*(n+1.0))
end

function get_all_coeffs_GS!()

    tab_YGS = Matrix{Float64}[]

    tab_lm = Any[]
    tab_index = zeros(Int64, lmax+1, mmax+1)


    println("Loading Gram Schmidt coefficients")

    index = 1
    nb_lm = 0

    for l=0:lmax
        for m=0:min(l,mmax)
            

            push!(tab_lm,[l,m])

            tab_index[l+1,m+1] = index

            index += 1

            nb_lm += 1

        end
    end

    

    for p=1:nb_lm

        l, m = tab_lm[p]

        n0 = 0
        if (m!=0)
            n0 = 1
        end

        namefile = "../../nb/GS_coeffs/h_1/l_"*string(l)*"/m_"*string(abs(m))*"/Y_l_"*string(l)*"_m_"*string(m)*".hdf5"
        file = h5open(namefile)
        data = read(file, "tabY")
        tab = [eval(Meta.parse(data[i-n0+1,j-n0+1])) for j=n0:nmax, i=n0:nmax]*1.0
        # tab[:,i+1] is the list of GS coefficients for the radial basis F^i
        # Flmn_GS(r) = sum_j={1}^{nmax-n0+1} tab[j,i+1] Flmj(r)

        push!(tab_YGS,tab)

        close(file)
    end

    return tab_YGS, tab_lm, tab_index
end

const tab_YGS, tab_lm, tab_index = get_all_coeffs_GS!()

function Flmn(l::Int64, m::Int64, n::Int64, xi::Float64)

    id = tab_index[l+1,abs(m)+1]
    sum = 0.0

    n0 = 0
    if (m != 0)
        n0 =1
    end

    for k=n0:nmax

        flmn = Flmn_bare(xi,l,m,k)
        sum += tab_YGS[id][k-n0+1,n-n0+1]*flmn 
    end

    return sum*sqrt(Delta)
end

function Flmn_x(l::Int64, m::Int64, n::Int64, x::Float64)

    if (x < 1.0)
        xi = (1.0+x)/(1.0-x)

        return Flmn(l,m,n,xi)
    else
        return 0.0
    end

end


function Dlmn_bare(xi::Float64, l::Int64, m::Int64, n::Int64)

    C0 = C0n(n)
    C1 = C1n(n)
    C2 = C2n(n)

    p = 1
    h = 1.0

    if (n > 1)

        pref1 = (xi-1)^(n-1)/(xi+h)^(n+p+2)
        term1 = C0 + C1*(xi-1) + C2*(xi-1)^2  

        pref2 = (xi-1)^(n-1)/(xi+h)^(p+n)
        term2 = l*(l+1)*(xi-1) + m^2/(xi+1)

        return (pref1*term1 - pref2*term2)
    elseif (n == 1)

        pref1 = 1.0/(xi+h)^(n+p+2)
        term1 = C0 + C1*(xi-1) + C2*(xi-1)^2  

        pref2 = 1.0/(xi+h)^(p+n)
        term2 = l*(l+1)*(xi-1) + m^2/(xi+1)

        return (pref1*term1 - pref2*term2)
    else # n = 0

        pref1 = 1.0/(xi+h)^(n+p+2)
        term1 =  C1 + C2*(xi-1)  

        pref2 = 1.0/(xi+h)^(p+n)
        term2 = l*(l+1) 

        return (pref1*term1 - pref2*term2)


    end
end

function Dlmn(l::Int64, m::Int64, n::Int64, xi::Float64)

    id = tab_index[l+1,abs(m)+1]
    sum = 0.0

    n0 = 0
    if (m != 0)
        n0 =1
    end

    for k=n0:nmax

        dlmn = Dlmn_bare(xi,l,m,k)
        sum += tab_YGS[id][k-n0+1,n-n0+1]*dlmn 
    end

    return sum*sqrt(Delta)
end

function Dlmn_x(l::Int64, m::Int64, n::Int64, x::Float64)

    
    if (x < 1.0)
        xi = (1.0+x)/(1.0-x)

        return Dlmn(l,m,n,xi)
    else
        return 0.0
    end

end


# Scalar product test (on x variable)


function scalar_productY(l::Int64, m::Int64, n::Int64, k::Int64, nbx::Int64=100)

    sum = 0.0

    for i=1:nbx 
        x = 1.0/nbx*(i-0.5)
        Fx =  Flmn_x(l,m,n,x)
        Dx =  Dlmn_x(l,m,k,x)

        sum += 2/(1-x)^2*Fx*Dx 
    end
    
    sum *= -(1.0/nbx)/Delta

    return sum 
end




# U(r,r') = -G/|r-r'|
# (r-r')^2 = (x-x')^2 + (y-y')^2 + (z-z')^2
# x = R cos phi 
# y = R sin phi 
# R = Delta sinh u sin v 
# z = Delta cosh u cos v
function interaction_potential_GS(u,v,phi,up,vp,phip)

    xi = cosh(u)
    xip = cosh(up)

    R = Delta * sinh(u) * sin(v)
    z = Delta * cosh(u) * cos(v)
    Rp = Delta * sinh(up) * sin(vp)
    zp = Delta * cosh(up) * cos(vp)

    x = R * cos(phi)
    y = R * sin(phi)
    xp = Rp * cos(phip)
    yp = Rp * sin(phip)

    diff_r = sqrt((x-xp)^2 + (y-yp)^2 + (z-zp)^2)

    println("U(r,r') [th.] = ",-G/diff_r)

    sum = 0.0


    nblm = length(tab_lm)

    for lm=1:nblm
        l, m = tab_lm[lm]

        n0 = 0
        if (m != 0)
            n0 = 1
        end

        for n=n0:nmax 



            # m >= 0
            psip_r = sqrt(4*pi*G)/Delta*Flmn(l,m,n,xi)*Ylm(l,m,v,phi)
            psip_rp = sqrt(4*pi*G)/Delta*Flmn(l,m,n,xip)*Ylm(l,m,vp,phip)

            sum += conj(psip_r)*psip_rp

            if (m != 0)

                psip_r = sqrt(4*pi*G)/Delta*Flmn(l,-m,n,xi)*Ylm(l,-m,v,phi)
                psip_rp = sqrt(4*pi*G)/Delta*Flmn(l,-m,n,xip)*Ylm(l,-m,vp,phip)
        
                sum += conj(psip_r)*psip_rp

                # println("sum = ",sum)
        
            end
        end
    end
    

    sum *= -1.0


    println("U(r,r') [SCF] = ",sum)

end


