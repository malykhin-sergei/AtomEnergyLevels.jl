using AtomEnergyLevels
using Test, Printf

function test01()
  @info "3D isotropic harmonic oscillator"
  @info "https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator"
  @info "Hamiltonian is:  H = -1/2μ⋅Δ + 1/2⋅k⋅r²"
  @info "Eigenvalues are: E = ħω⋅(2nᵣ + l + 3/2),"
  @info "where nᵣ = 0, 1, 2...; l = 0, 1, 2..., and ω² = k/μ"
  r, n, dx = begin
      x = range(-30, 20, length = 501)
      exp.(x), x.len, x.step.hi
  end
  conf = ((0,0,0,), (0,0,0), (0,0,0))
  ψ = radial_shr_eq((r, n, dx), 1/2*r.^2, conf).orbitals
  @info "nᵣ\tl\tϵ(calc.)\tϵ(exact)\tΔϵ"
  nlevels = length(collect(Iterators.flatten(conf)))
  ϵ_exact = zeros(nlevels)
  k = 1
  for (i, subshell) in enumerate(conf)
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
  return abs_error
end

function test02()
  @info "Hooke's atom"
  @info "https://en.wikipedia.org/wiki/Hooke's_atom"
  @info "An exact solution for Hooke's atom is E = 2.0 a.u."
  @info "For Slater Xα method α is adjustable parameter."
  @info "Here we reproduce exact Hooke atom energy with"
  @info "α = 0.83685294"
  r, n, dx = begin
      x = range(-30, 20, length = 301)
      exp.(x), x.len, x.step.hi
  end
  E = lda((r, n, dx), conf = 2, xc = x -> LDA_X(x; α = 0.83685294), vp = 1/8 * r.^2, β = 0.8).energy
  return E-2.0
end

function test03()
  @info "Argon atom"
  @info "https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-16"
  @info "Etot = -525.946195 (NIST)"
  E, ρ, ψ = lda(conf = atomic_electron_configuration[:Ar], β=0.6)
  @info "Energy levels for the Ar atom are:"
  for (i, ϵᵢ) in enumerate(ψ[3,:])
    l = ψ[1, i]
    n = l + ψ[2, i] + 1
    @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f\n", n, atomic_shell[l+1], ψ[4, i], ϵᵢ)
  end
  return -525.946195-E
end

@testset "AtomEnergyLevels.jl" begin
  @test test01() ≈ 0.0 atol = 1e-10
  @test test02() ≈ 0.0 atol = 5e-7
  @test test03() ≈ 0.0 atol = 5e-7
end
