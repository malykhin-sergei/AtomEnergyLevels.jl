using AtomEnergyLevels
using Test, Printf

@testset "AtomEnergyLevels.jl" begin
  @testset "Differentiation matrices" begin
    x = -6:0.01:6
    d²f_ex = 2(2x.^2 .- 1).*exp.(-x.^2)
    d²f_ap = laplacian(x) * exp.(-x.^2)
    @test d²f_ex ≈ d²f_ap
  end

  @testset "Parsing electronic configuration" begin
    @test c"1s1"              == ((1.0,),)
    @test c"[He] 2s1"         == ((2.0, 1.0),)
    @test c"[He] 2s2 2p4"     == ((2.0, 2.0), (4.0,))
    @test c"[Xe] 4f1 5d1 6s2" == ((2.0, 2.0, 2.0, 2.0, 2.0, 2.0), (6.0, 6.0, 6.0, 6.0), (10.0, 10.0, 1.0), (1.0,))
    @test_throws DomainError conf_enc("[Ge] 1s1")
    @test_throws DomainError conf_enc("[Rn]", maxn = 5)
    @test_throws DomainError conf_enc("[He] 2d1")
    @test_throws DomainError conf_enc("1s2", spins = 1)
    @test_throws DomainError conf_enc("[He] 1s2")
  end
  
  @testset "DFT functionals" begin
    rₛ = [2, 5, 10, 20, 50, 100]
    ρ = 1 ./ (4π/3 * rₛ .^ 3)

    ϵ_vwn5 = [-0.0447827886146218,  -0.0281337622897314,  -0.0185445271694030,  
              -0.0115476822648961,  -0.0057034884827245,  -0.0031846468815323]
    v_vwn5 = [-0.0516038239497904,  -0.0333841710353646,  -0.0225183261458631,  
              -0.0143300153523442,  -0.0072495456200668,  -0.0041038158928284]

    ϵ_chac = [-0.0434298249201931,  -0.0276172696449811,  -0.0183235042675089,  
              -0.0113396589743083,  -0.0054215811399318,  -0.0029196229936562]
    v_chac = [-0.0499160959362279,  -0.0326396845509736,  -0.0222371196290851,  
              -0.0141507122981238,  -0.0069772050422708,  -0.0038156875497190]

    for (i, ρᵢ) in enumerate(ρ)
      vᵢ, ϵᵢ = LDA_C_VWN(ρᵢ)
      @test ϵᵢ ≈ ϵ_vwn5[i]
      @test vᵢ ≈ v_vwn5[i]
      vᵢ, ϵᵢ = LDA_C_CHACHIYO(ρᵢ)
      @test ϵᵢ ≈ ϵ_chac[i]
      @test vᵢ ≈ v_chac[i]
    end
  end

  @testset "3D isotropic harmonic oscillator" begin
    @info "3D isotropic harmonic oscillator"
    @info "https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator"
    @info "Hamiltonian is:  H = -1/2μ⋅Δ + 1/2⋅k⋅r²"
    @info "Eigenvalues are: E = ħω⋅(2nᵣ + l + 3/2),"
    @info "where nᵣ = 0, 1, 2...; l = 0, 1, 2..., and ω² = k/μ"
    config = ((0,0,0,), (0,0,0), (0,0,0))
    ψ = radial_shr_eq(r -> 1/2*r^2, conf = config).orbitals
    @info "nᵣ\tl\tϵ(calc.)\tϵ(exact)\tΔϵ"
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      ϵ_calc = orbital.ϵᵢ
      ϵ_exact = 2nᵣ + l + 3/2      
      @info @sprintf("%i\t%s\t%10.8f\t%10.8f\t%+0.6e",
      nᵣ, shells[l], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)
      @test ϵ_calc ≈ ϵ_exact atol = 1e-10
    end
  end

  @testset "Hooke's atom energy" begin
    @info "Hooke's atom"
    @info "https://en.wikipedia.org/wiki/Hooke's_atom"
    @info "An exact solution for Hooke's atom is E = 2.0 a.u."
    @info "For the Xα method α is adjustable parameter."
    @info "Here we reproduce exact Hooke atom energy with α = 0.83685294"
    Etot = lda(2, conf = c"[He]", 
               xc! = (ρ, vxc, exc) -> Xα!(ρ, vxc, exc, α = 0.83685294), 
               Vex = r -> 1/8 * r^2).energy.total
    @test Etot ≈ 2.0 atol = 1e-7
  end

  @testset "Argon atom" begin
    @info "Argon atom"
    @info "https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-16"
    @info "Etot = -525.946195 (NIST)"
    E, ρ, ψ = lda(18)
    @info "Energy levels for the Ar atom are:"
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      nᵢ, ϵᵢ = orbital.nᵢ, orbital.ϵᵢ
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, shells[l], nᵢ, ϵᵢ)
    end
    @test E.total     ≈  -525.946195 atol = 5e-7 # Etot
    @test E.kinetic   ≈   524.969812 atol = 1e-6 # Ekin
    @test E.hartree   ≈   231.458124 atol = 1e-6 # Ecoul
    @test E.potential ≈ -1253.131982 atol = 1e-6 # Eenuc 
    @test E.xc        ≈   -29.242149 atol = 1e-6 # Exc

    @test ψ[(nᵣ=0, l=0)].ϵᵢ ≈  -113.800134 atol = 1e-6 # 1s
    @test ψ[(nᵣ=1, l=0)].ϵᵢ ≈   -10.794172 atol = 1e-6 # 2s
    @test ψ[(nᵣ=0, l=1)].ϵᵢ ≈    -8.443439 atol = 1e-6 # 2p
    @test ψ[(nᵣ=2, l=0)].ϵᵢ ≈    -0.883384 atol = 1e-6 # 3s
    @test ψ[(nᵣ=1, l=1)].ϵᵢ ≈    -0.382330 atol = 1e-6 # 3p
  end
  
  @testset "Ar+" begin
    @info "Ar+ [Ne] 3s2 3p5"
    @info "https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-16"
    @info "Etot = -525.351708  (NIST)"
    E, ρ, ψ = lda(18, conf = c"[Ne] 3s2 3p5")
    @info "Energy levels for Ar+ are:"
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      nᵢ, ϵᵢ = orbital.nᵢ, orbital.ϵᵢ
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, shells[l], nᵢ, ϵᵢ)
    end
    @test E.total     ≈  -525.351708 atol = 5e-7 # Etot
    @test E.kinetic   ≈   524.405209 atol = 1e-6 # Ekin
    @test E.hartree   ≈   222.915201 atol = 1e-6 # Ecoul
    @test E.potential ≈ -1243.839462 atol = 1e-6 # Eenuc 
    @test E.xc        ≈   -28.832655 atol = 1e-6 # Exc
    
    @test ψ[(nᵣ=0, l=0)].ϵᵢ ≈  -114.320786 atol = 1e-6 # 1s
    @test ψ[(nᵣ=1, l=0)].ϵᵢ ≈   -11.303467 atol = 1e-6 # 2s
    @test ψ[(nᵣ=0, l=1)].ϵᵢ ≈    -8.954227 atol = 1e-6 # 2p
    @test ψ[(nᵣ=2, l=0)].ϵᵢ ≈    -1.337427 atol = 1e-6 # 3s
    @test ψ[(nᵣ=1, l=1)].ϵᵢ ≈    -0.816635 atol = 1e-6 # 3p
  end
end
