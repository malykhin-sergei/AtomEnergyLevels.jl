using AtomEnergyLevels
using Test, Printf

@testset "AtomEnergyLevels.jl" begin
  @testset "3D isotropic harmonic oscillator" begin
    @info "3D isotropic harmonic oscillator"
    @info "https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator"
    @info "Hamiltonian is:  H = -1/2μ⋅Δ + 1/2⋅k⋅r²"
    @info "Eigenvalues are: E = ħω⋅(2nᵣ + l + 3/2),"
    @info "where nᵣ = 0, 1, 2...; l = 0, 1, 2..., and ω² = k/μ"
    config = ((0,0,0,), (0,0,0), (0,0,0))
    ψ = radial_shr_eq(r -> 1/2*r^2, conf = config).orbitals
    @info "nᵣ\tl\tϵ(calc.)\tϵ(exact)\tΔϵ"
    nlevels = length(collect(Iterators.flatten(config)))
    ϵ_exact = zeros(nlevels)
    k = 1
    for (i, subshell) in enumerate(config)
      l = i - 1
      for (j, nᵢ) in enumerate(subshell)
        nᵣ = j - 1
        ϵ_exact[k] = 2nᵣ + l + 3/2
        @info @sprintf("%i\t%s\t%10.8f\t%10.8f\t%+0.6e",
        nᵣ, atomic_shell[i], ψ[3,k], ϵ_exact[k], ϵ_exact[k] - ψ[3,k])
        k += 1
      end
    end
    abs_error = maximum(abs.(ϵ_exact .- ψ[3,:]))
    @info @sprintf("Maximum absolute error is: max|Δϵ| = %0.6e", abs_error)
    @test abs_error ≈ 0.0 atol = 1e-10
  end

  @testset "Hooke's atom energy" begin
    @info "Hooke's atom"
    @info "https://en.wikipedia.org/wiki/Hooke's_atom"
    @info "An exact solution for Hooke's atom is E = 2.0 a.u."
    @info "For the Xα method α is adjustable parameter."
    @info "Here we reproduce exact Hooke atom energy with α = 0.83685294"
    Etot = lda(2, conf = 2, xc = x -> LDA_X(x; α = 0.83685294), 
                Vex = r -> 1/8 * r^2, β = 0.8).energy.total
    @test Etot ≈ 2.0 atol = 1e-7
  end

  @testset "Argon atom" begin
    @info "Argon atom"
    @info "https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-16"
    @info "Etot = -525.946195 (NIST)"
    E, ρ, ψ = lda(18, β=0.6)
    @info "Energy levels for the Ar atom are:"
    p = sortperm(ψ[3,:])
    ψ = ψ[:, p]
    for (i, ϵᵢ) in enumerate(ψ[3,:])
      l, nᵣ, nᵢ = ψ[1, i], ψ[2, i], ψ[4, i]
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, atomic_shell[l+1], nᵢ, ϵᵢ)
    end
    @test E.total     ≈  -525.946195 atol = 5e-7 # Etot
    @test E.kinetic   ≈   524.969812 atol = 1e-6 # Ekin
    @test E.hartree   ≈   231.458124 atol = 1e-6 # Ecoul
    @test E.potential ≈ -1253.131982 atol = 1e-6 # Eenuc 
    @test E.xc        ≈   -29.242149 atol = 1e-6 # Exc
    @test ψ[3, 1]     ≈  -113.800134 atol = 1e-6 # 1s
    @test ψ[3, 2]     ≈   -10.794172 atol = 1e-6 # 2s
    @test ψ[3, 3]     ≈    -8.443439 atol = 1e-6 # 2p
    @test ψ[3, 4]     ≈    -0.883384 atol = 1e-6 # 3s
    @test ψ[3, 5]     ≈    -0.382330 atol = 1e-6 # 3p
  end
  
  @testset "Ar+" begin
    @info "Ar+ [Ne] 3s2 3p5"
    @info "https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-16"
    @info "Etot = -525.351708  (NIST)"
    E, ρ, ψ = lda(18, conf = conf_enc("[Ne] 3s2 3p5"), β=0.6)
    @info "Energy levels for Ar+ are:"
    p = sortperm(ψ[3,:])
    ψ = ψ[:, p]
    for (i, ϵᵢ) in enumerate(ψ[3,:])
      l, nᵣ, nᵢ = ψ[1, i], ψ[2, i], ψ[4, i]
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, atomic_shell[l+1], nᵢ, ϵᵢ)
    end
    @test E.total     ≈  -525.351708 atol = 5e-7 # Etot
    @test E.kinetic   ≈   524.405209 atol = 1e-6 # Ekin
    @test E.hartree   ≈   222.915201 atol = 1e-6 # Ecoul
    @test E.potential ≈ -1243.839462 atol = 1e-6 # Eenuc 
    @test E.xc        ≈   -28.832655 atol = 1e-6 # Exc
    @test ψ[3, 1]     ≈  -114.320786 atol = 1e-6 # 1s
    @test ψ[3, 2]     ≈   -11.303467 atol = 1e-6 # 2s
    @test ψ[3, 3]     ≈    -8.954227 atol = 1e-6 # 2p
    @test ψ[3, 4]     ≈    -1.337427 atol = 1e-6 # 3s
    @test ψ[3, 5]     ≈    -0.816635 atol = 1e-6 # 3p
  end
end
