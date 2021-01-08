using AtomEnergyLevels, Libxc

function SLAPW92!(ρ, vxc, εxc)
    lda_x = Functional(:lda_x)
    lda_c = Functional(:lda_c_pw)
    exch = evaluate(lda_x, rho=ρ)
    corr = evaluate(lda_c, rho=ρ)
    for i in eachindex(ρ)
        vxc[i] = exch.vrho[i] + corr.vrho[i]
        εxc[i] = exch.zk[i] + corr.zk[i]
    end
end

function atom_ground_state(z, xcfun)
    ψ = lda(z, xc! = xcfun).orbitals
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      nᵢ, ϵᵢ = orbital.nᵢ, orbital.ϵᵢ
      n = nᵣ + l + 1
      @printf("\t%i%s\t(%4.1f)\t%14.6f\n", n, shells[l], nᵢ, ϵᵢ)
    end
end
  
atom_ground_state(18, SLAPW92!)
