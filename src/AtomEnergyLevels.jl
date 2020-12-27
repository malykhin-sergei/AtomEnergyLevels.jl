module AtomEnergyLevels

include("integrate.jl")
include("sincdif.jl")
include("conf_parse.jl")
include("periodic_table.jl")
include("dft_xc_functionals.jl")
include("rshreq.jl")
include("TF.jl")

import LinearAlgebra: eigen!, diagind, Diagonal, Symmetric
import Printf: @sprintf

export laplacian, radial_shr_eq, TF, lda
export LDA_X, LDA_C_CHACHIYO, LDA_C_VWN
export SVWN!, Xα!
export conf_enc, Nₑ
export atomic_shell, atomic_electron_configuration
export laplacian
export radial_shr_eq

"""
```julia
function lda(Z, x = -30.0:0.1:20.0; conf = atomic_electron_configuration[Z],
             xc! = SVWN!, Vex = r -> -Z / r,
             ρᵢₙ = TF.(exp.(x), Nₑ(conf)) .* exp.(x), 
             β = 0.3, δn = 1.0e-6, δE = 5.0e-7, maxit = 100,
             μ = 1, α = 1e5)
```
Solve Kohn-Sham DFT Self-Consistent Field equations for an atom using
[local density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

On input:

* `Z`     - nuclear charge
* `x`     - grid, rᵢ = exp(xᵢ)
* `conf`  - electronic configuration
* `xc`    - LDA exchange-correlation functional
* `Vex`   - external potential
* `ρᵢₙ`   - input density (if nothing provided, Thomas-Fermi density is used)
* `β`     - density mixing parameter, ρₙ = (1 - β) × ρₙ₋₁ + β × ρₙ
* `δn`    - density convergence criterion
* `δE`    - energy convergence criterion
* `maxit` - maximum permissible number of SCF iterations

On output:

*  `lda(...).energy.total` - total energy, a.u.
*  `lda(...).density`  - electron density ρ(x) = 1/4π * ∑nᵢ*ψᵢ²(x)
*  `lda(...).orbitals` - matrix, containing complete set of orbitals ψᵢ(x), corresponding eigenvalues εᵢ (energy levels), azimuthal and radial quantum numbers l, nᵣ, and level populations nᵢ (as listed in `conf`).

# Example
[Hooke atom](https://en.wikipedia.org/wiki/Hooke's_atom)

```julia-repl
julia> lda(2, Vex = r -> 1/8 * r^2, β = 0.8).energy.total;
[ Info: Using logarithmic 501 point grid with step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle           energy          |Δρ|
[ Info:   0           2.103849      1.895566
[ Info:   1           2.013080      0.357780
[ Info:   2           2.022087      0.073897
[ Info:   3           2.025397      0.014779
[ Info:   4           2.026076      0.003000
[ Info:   5           2.026201      0.000621
[ Info:   6           2.026224      0.000130
[ Info:   7           2.026228      0.000028
[ Info:   8           2.026229      0.000006
[ Info:   9           2.026229      0.000001
[ Info:  10           2.026229      0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022579
│       EXCHANGE-CORRELATION          -0.523773
│       ELECTRON-NUCLEAR               0.899965
│       TOTAL ENERGY                   2.026229
└       VIRIAL RATIO                  -2.229260

```
Compare with exact Kohn-Sham values from Table II of S. Kais et al
// Density functionals and dimensional renormalization for an exactly solvable model,
JCP 99, 417 (1993); `http://dx.doi.org/10.1063/1.465765`
```
        ELECTRON KINETIC               0.6352
        ELECTRON-ELECTRON              1.0320
        EXCHANGE-CORRELATION          -0.5553
        ELECTRON-NUCLEAR               0.8881
        TOTAL ENERGY                   2.0000
        VIRIAL RATIO                  -2.1486
```
"""
function lda(Z,
             x = -30.0:0.1:20.0;             
          conf = atomic_electron_configuration[Z],
           xc! = SVWN!,
           Vex = r -> -Z / r,
           ρᵢₙ = TF.(exp.(x), Nₑ(conf)) .* exp.(x),
             β = 0.3,
            δn = 1.0e-6,
            δE = 5.0e-7,
         maxit = 100,
             μ = 1,
             α = 1e5)

  # radial grid
  r, n, dx = exp.(x), length(x), step(x)
  r², sqr = exp.(2x), exp.(x/2)

  @info @sprintf("Using logarithmic %3i point grid with step dx = %5.4f", n, dx)

  vp = Vex.(r)

  Δ = laplacian(x)
  L = Δ - Diagonal(fill(1/4, n))

  Eₜₒₜ = 0.0; Δρ = 0.0; Q = Nₑ(conf)
  ρₒᵤₜ, V, vh, vxc, εxc = [similar(ρᵢₙ) for _ = 1:5]

  @info @sprintf("Starting SCF procedure with density mixing parameter β = %5.4f
      and convergence threshold |Δρ| ≤ %e", β, δn)
  @info "cycle\t\tenergy\t\t|Δρ|"
  for i = 0:maxit
    # Solve the Poisson equation to obtain Hartree potential
    vh = L \ (-4π * ρᵢₙ .* sqr .* r)

    # Apply boundary conditions at r → 0 and r → ∞
    vh .-= (vh[n] - Q / sqr[n])/sqr[n] .* sqr .+
           (vh[1] - Q * sqr[1])*sqr[1] ./ sqr
    # Change variable vh(x) → vh(r)
    vh ./= sqr

    # Calculate exchange-correlation potential and energy density
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#the-xc-term)
    xc!(ρᵢₙ ./ r, vxc, εxc)

    # Assemble Kohn-Sham potential
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#kohn-sham-equations)
    V = vp + vh + vxc

    # Solve the Schrödinger equation to find new density
    # and bands energy ∑nᵢεᵢ
    H = -1/2μ*Δ + Diagonal(V .* r²)

    ∑nᵢεᵢ = 0.0; ρₒᵤₜ .= 0.0

    # the equation is solved separately for each subshell: s, p, d, f
    for (l, subshell) in enumerate(conf)
      Hl = H + Diagonal(fill(1/2μ * (l - 1/2)^2, n))
      θ, y = eigen!(Symmetric(Hl), Symmetric(Hl + α*Diagonal(r²)))
      ε = α*θ ./ (1 .- θ)
      for (nᵣ, nᵢ) in enumerate(subshell)
        y[:, nᵣ] /= sqrt(∫(dx, y[:, nᵣ] .^ 2 .* r²))
        ρₒᵤₜ .+= nᵢ / 4π * y[:, nᵣ] .^ 2
        ∑nᵢεᵢ += nᵢ * ε[nᵣ]
      end
    end

    E = Eₜₒₜ

    # DFT total energy:
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#total-energy)
    Eₜₒₜ = ∑nᵢεᵢ + 4π * ∫(dx, ρₒᵤₜ .* (-1/2 * vh .- vxc .+ εxc) .* r²)

    Δρ = 4π * ∫(dx, abs.(ρₒᵤₜ - ρᵢₙ) .* r²)
    @info @sprintf "%3i\t%14.6f\t%12.6f\n" i Eₜₒₜ Δρ

    # density converged if there are no more charge oscillations
    Δρ < δn && abs(Eₜₒₜ - E) < δE && break

    # admix density from previous iteration to converge SCF
    ρᵢₙ = (1.0 - β) * ρᵢₙ + β * ρₒᵤₜ
  end
  Δρ > δn && @warn "SCF does not converged after $maxit iterations!"

  # obtain orbitals
  ∑ε, ρ, ψ = radial_shr_eq(V, x, conf = conf, μ = μ)
  
  # remove numerical noise
  map!(x -> ifelse(x > eps(), x, 0.0), ρ, ρ)

  # total energy components
  Ek  = ∑ε - 4π * ∫(dx, ρ .* V .* r²)
  Eh  = 2π * ∫(dx, ρ .* vh .* r²)
  Exc = 4π * ∫(dx, ρ .* εxc .* r²)
  Ep  = 4π * ∫(dx, ρ .* vp .* r²)
  E   = Ek + Ep + Eh + Exc
  VR  = (Eh + Exc + Ep) / -Ek

  @info @sprintf("RESULTS SUMMARY:
      ELECTRON KINETIC       %16.6f
      ELECTRON-ELECTRON      %16.6f
      EXCHANGE-CORRELATION   %16.6f
      ELECTRON-NUCLEAR       %16.6f
      TOTAL ENERGY           %16.6f
      VIRIAL RATIO           %16.6f\n", Ek, Eh, Exc, Ep, E, VR)
      
  return (energy = (total = E, kinetic = Ek, hartree = Eh, 
          xc = Exc, potential = Ep), density = ρ, orbitals = ψ)
end
end # module
