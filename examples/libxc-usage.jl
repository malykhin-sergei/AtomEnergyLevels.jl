using AtomEnergyLevels, Libxc, Printf, Test

function libxcfun(xc::Symbol, threshold = eps())
    fxc = Functional(xc)
    if getproperty(fxc, :kind) ≠ :exchange_correlation
        throw(DomainError(getproperty(fxc, :identifier), "is not an exchange-correlation functional"))
    end
    fxc.density_threshold = threshold
    function xc!(ρ, vxc, exc)
        xc = evaluate(fxc, rho = ρ)
        for i in eachindex(ρ)
            @inbounds vxc[i], exc[i] = xc.vrho[i], xc.zk[i]
        end
        return vxc, exc
    end
end

function libxcfun(x::Symbol, c::Symbol, threshold = eps())
    fx, fc = Functional(x), Functional(c)
    fx.density_threshold = fc.density_threshold = threshold
    function xc!(ρ, vxc, exc)
        exch = evaluate(fx, rho = ρ)
        corr = evaluate(fc, rho = ρ)
        for i in eachindex(ρ)
            vxc[i] = exch.vrho[i] + corr.vrho[i]
            exc[i] = exch.zk[i]   + corr.zk[i]
        end
        return vxc, exc
    end
end 

function atom_ground_state(z, xc_func)
    results = lda(z, xc_func! = libxcfun(xc_func...))
    ψ = results.orbitals
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      nᵢ, ϵᵢ = orbital.nᵢ, orbital.ϵᵢ
      n = nᵣ + l + 1
      @printf("\t%i%s\t(%4.1f)\t%14.6f\n", n, shells[l + 1], nᵢ, ϵᵢ)
    end
    return results.energy.total
end
  
atom_ground_state(18, (:lda_x, :lda_c_vwn))
atom_ground_state(18, (:lda_x, :lda_c_pz))
atom_ground_state(18, (:lda_x, :lda_c_pw))
atom_ground_state(18, (:lda_xc_ksdt,))
atom_ground_state(18, (:lda_xc_gdsmfb,))
atom_ground_state(18, (:lda_xc_teter93,))