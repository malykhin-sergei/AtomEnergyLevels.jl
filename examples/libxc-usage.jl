using AtomEnergyLevels, Libxc, Printf

function libxcfun(xc::Symbol)
    fxc = Functional(xc)

    if getproperty(fxc, :kind) ≠ :exchange_correlation
        throw(DomainError(getproperty(fxc, :identifier), "is not an exchange-correlation functional"))
    end

    function xc_lda(ρ, vxc = similar(ρ), exc = similar(ρ))
        if !(length(ρ) == length(vxc) == length(exc))
            throw(DimensionMismatch("Dimensions of the arrays ρ, vxc, and exc doesn't match"))
        end
        evaluate!(fxc, rho = ρ, vrho = vxc, zk = exc)
        return vxc, exc
    end

    if getproperty(fxc, :family) == :lda
        return XC_Functional(xc_type_lda, xc_lda)
    else
        return nothing
    end
end

function libxcfun(x::Symbol, c::Symbol)
    fx, fc = Functional(x), Functional(c)

    if getproperty(fx, :kind) ≠ :exchange || getproperty(fc, :kind) ≠ :correlation
        throw(DomainError(getproperty(fxc, :identifier), "is not an exchange-correlation functional"))
    end

    function xc_lda(ρ, vxc = similar(ρ), exc = similar(ρ))
        if !(length(ρ) == length(vxc) == length(exc))
            throw(DimensionMismatch("Dimensions of the arrays ρ, vxc, and exc doesn't match"))
        end
        exch = evaluate(fx, rho = ρ)
        corr = evaluate(fc, rho = ρ)
        @inbounds for i in eachindex(ρ)
            vxc[i] = exch.vrho[i] + corr.vrho[i]
            exc[i] = exch.zk[i]   + corr.zk[i]
        end
        return vxc, exc  
    end

    if getproperty(fx, :family) == :lda && getproperty(fc, :family) == :lda
        return XC_Functional(xc_type_lda, xc_lda)
    else
        return nothing
    end
end

lda(2, xc = libxcfun(:lda_x, :lda_c_vwn))

function zoo(z)
    function is_lda_xc(x::Functional)
        id = getproperty(x, :identifier)
        occursin(r"_\dd_", string(id)) && return false

        getproperty(x, :family) == :lda &&
        getproperty(x, :kind)   == :exchange_correlation &&
    
        intersect([:vxc, :exc], getproperty(x, :flags)) == [:vxc, :exc]
    end

    function is_lda_c(x::Functional)
        id = getproperty(x, :identifier)
        occursin(r"_\dd_", string(id)) && return false

        getproperty(x, :family) == :lda &&
        getproperty(x, :kind)   == :correlation &&
    
        intersect([:vxc, :exc], getproperty(x, :flags)) == [:vxc, :exc]
    end

    xc_funcs = [] # filter(is_lda_xc, Functional.(available_functionals()))
     c_funcs = filter(is_lda_c, Functional.(available_functionals()))

    ids_xc = map(λ -> (getproperty(λ, :identifier),), xc_funcs)
    ids_c  = map(λ -> (:lda_x, getproperty(λ, :identifier)), c_funcs)

    energies = []

    for xc in vcat(ids_xc, ids_c)
        @printf("Testing %s\n", string(xc))
        try
            push!(energies, lda(z, xc = libxcfun(xc...)).energy.total)
        catch
            push!(energies, NaN)
            continue
        end
    end

    @printf("== Total energies of the %s atom in various LDA functionals ==\n", keys(atom)[z])
    for (xc, en) in zip(vcat(xc_funcs, c_funcs), energies)
        if getproperty(xc, :kind) == :correlation
            abbr = string.("lda_x + ", getproperty(xc, :identifier))
        else
            abbr = string.("        ", getproperty(xc, :identifier))
        end
        name = getproperty(xc, :name)
        year = last(split(getproperty(xc, :references)[1].reference))
        @printf("%-30s\t%12.6f\t%s %s\n", abbr, en, name, year)
    end
end

zoo(2)
#=
== Total energies of the He atom in various LDA functionals ==
        lda_xc_teter93             -2.833976    Teter 93 (1996)
        lda_xc_zlp                 -2.913789    Zhao, Levy & Parr, Eq. (20) (1993)
        lda_xc_ksdt                -2.833100    Karasiev, Sjostrom, Dufty & Trickey (2014)
        lda_xc_lp_a                -2.869133    Lee-Parr reparametrization A (1990)
        lda_xc_lp_b                -2.903489    Lee-Parr reparametrization B (1990)
        lda_xc_gdsmfb              -2.832794    Groth, Dornheim, Sjostrom, Malone, Foulkes, Bonitz (2017)
        lda_xc_bn05                -2.244445    Baer and Neuhauser, gamma=1 (2005)
lda_x + lda_c_wigner               -2.818849    Wigner (1938)
lda_x + lda_c_rpa                  -2.837328    Random Phase Approximation (RPA) (1957)
lda_x + lda_c_hl                   -2.839923    Hedin & Lundqvist (1971)
lda_x + lda_c_gl                   -2.860137    Gunnarson & Lundqvist (1976)
lda_x + lda_c_xalpha               -3.170112    Slater's Xalpha (1951)
lda_x + lda_c_vwn                  -2.834836    Vosko, Wilk & Nusair (VWN5) (1980)
lda_x + lda_c_vwn_rpa              -2.872169    Vosko, Wilk & Nusair (VWN5_RPA) (1980)
lda_x + lda_c_pz                   -2.834291    Perdew & Zunger (1981)
lda_x + lda_c_pz_mod               -2.834313    Perdew & Zunger (Modified) parts
lda_x + lda_c_ob_pz                -2.830831    Ortiz & Ballone (PZ parametrization) (1994)
lda_x + lda_c_pw                   -2.834455    Perdew & Wang (1992)
lda_x + lda_c_pw_mod               -2.834455    Perdew & Wang (modified) (http://dft.rutgers.edu/pubs/PBE.asc)
lda_x + lda_c_ob_pw                -2.831420    Ortiz & Ballone (PW parametrization) (1994)
lda_x + lda_c_vbh                  -2.870439    von Barth & Hedin (1972)
lda_x + lda_c_ml1                  -2.779743    Modified LSD (version 1) of Proynov and Salahub (1994)
lda_x + lda_c_ml2                  -2.751033    Modified LSD (version 2) of Proynov and Salahub (1994)
lda_x + lda_c_gombas               -2.837973    Gombas (1965)
lda_x + lda_c_pw_rpa               -2.871024    Perdew & Wang (fit to the RPA energy) (1992)
lda_x + lda_c_rc04                 -2.816964    Ragot-Cortona (2004)
lda_x + lda_c_vwn_1                -2.834836    Vosko, Wilk & Nusair (VWN1) (1980)
lda_x + lda_c_vwn_2                -2.834836    Vosko, Wilk & Nusair (VWN2) (1980)
lda_x + lda_c_vwn_3                -2.834836    Vosko, Wilk & Nusair (VWN3) (1980)
lda_x + lda_c_vwn_4                -2.834836    Vosko, Wilk & Nusair (VWN4) (1980)
lda_x + lda_c_chachiyo             -2.831427    Chachiyo simple 2 parameter correlation (2016)
lda_x + lda_c_lp96                       NaN    Liu-Parr correlation (1996)
lda_x + lda_c_chachiyo_mod         -2.831427    Chachiyo simple 2 parameter correlation with modified spin scaling [cond-mat.mtrl-sci]
lda_x + lda_c_karasiev_mod         -2.833235    Karasiev reparameterization of Chachiyo [cond-mat.mtrl-sci]
lda_x + lda_c_mcweeny              -2.850426    McWeeny 76 3--31
lda_x + lda_c_br78                 -2.764565    Brual & Rothstein 78 (1978)
lda_x + lda_c_pk09                 -2.834520    Proynov and Kong 2009 (2009)
lda_x + lda_c_ow_lyp               -2.780839    Wigner with corresponding LYP parameters (1995)
lda_x + lda_c_ow                   -2.784822    Optimized Wigner (1995)
lda_x + lda_c_gk72                 -2.839369    Gordon and Kim 1972 (1972)
lda_x + lda_c_karasiev             -2.833235    Karasiev reparameterization of Chachiyo (2016)
lda_x + lda_c_upw92                -2.834749    Ruggeri, Rios, and Alavi unrestricted fit (2018)
lda_x + lda_c_rpw92                -2.834455    Ruggeri, Rios, and Alavi restricted fit (2018)
=#
