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
             xc! = SVWN!, Vex = r -> -Z / r, δn = 1e-7, maxit = 100, μ = 1)
```
Solve Kohn-Sham DFT Self-Consistent Field equations for an atom using
[local density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

On input:

* `Z`     - nuclear charge
* `x`     - grid, rᵢ = exp(xᵢ)
* `conf`  - electronic configuration
* `xc`    - LDA exchange-correlation functional
* `Vex`   - external potential
* `δn`    - density convergence criterion
* `maxit` - maximum permissible number of SCF iterations

On output:

*  `lda(...).energy.total` - total energy, a.u.
*  `lda(...).density`  - electron density ρ(x) = 1/4π * ∑nᵢ*ψᵢ²(x)
*  `lda(...).orbitals` - orbitals ψᵢ(x), energy levels εᵢ, azimuthal and radial quantum numbers l, nᵣ, and populations nᵢ (as listed in `conf`).

# Example
[Hooke atom](https://en.wikipedia.org/wiki/Hooke's_atom)

```julia-repl
julia> lda(2, Vex = r -> 1/8 * r^2).energy.total;
[ Info: Using logarithmic 501 point grid with step dx = 0.1000
[ Info: Starting SCF procedure with convergence threshold |Δρ| ≤ 1.000000e-07
[ Info: cycle           energy          |Δρ|
[ Info:   0           2.103849      1.895566
[ Info:   1           2.013080      0.357780
[ Info:   2           2.023521      0.039644
[ Info:   3           2.026081      0.003607
[ Info:   4           2.026212      0.000404
[ Info:   5           2.026228      0.000045
[ Info:   6           2.026229      0.000005
[ Info:   7           2.026229      0.000001
[ Info:   8           2.026229      0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022579
│       EXCHANGE-CORRELATION          -0.523773
│       ELECTRON-NUCLEAR               0.899965
│       TOTAL ENERGY                   2.026229
└       VIRIAL RATIO                  -2.229260
```
"""
function lda(Z,
             x = -30.0:0.1:20.0;             
          conf = atomic_electron_configuration[Z],
           xc! = SVWN!,
           Vex = r -> -Z / r,
            δn = 1e-7,
         maxit = 100,
             μ = 1)
  # 
  β = 0.8; α = 1e5

  # radial grid
  r, n, dx = exp.(x), length(x), step(x)
  r², sqr = exp.(2x), exp.(x/2)

  @info @sprintf("Using logarithmic %3i point grid with step dx = %5.4f", n, dx)

  Q = Nₑ(conf)
  ρᵢₙ = TF.(r, Q) .* r
  vp = Vex.(r)

  Δ = laplacian(x)
  L = Δ - Diagonal(fill(1/4, n))

  Eₜₒₜ = 0.0; Δρ = 0.0; 
  
  ρₒᵤₜ, prev_ρᵢₙ, prev_ρₒᵤₜ, V, vh, vxc, εxc = [similar(ρᵢₙ) for _ = 1:7]

  @info @sprintf("Starting SCF procedure with convergence threshold |Δρ| ≤ %e", δn)
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
    @info @sprintf "%3i\t%14.6f\t%12.6f" i Eₜₒₜ Δρ

    # density converged if there are no more charge oscillations
    Δρ < δn && break

    if i > 0
      dρₒᵤₜ = ∫(dx, (ρₒᵤₜ - prev_ρₒᵤₜ) .* r) / ∫(dx, (ρᵢₙ - prev_ρᵢₙ) .* r)
      β = min(max(1/(1 - dρₒᵤₜ), 0.1), 0.9)
    end

    prev_ρᵢₙ  .= ρᵢₙ
    prev_ρₒᵤₜ .= ρₒᵤₜ

    # admix density from previous iteration to converge SCF
    ρᵢₙ = (1 - β) * ρᵢₙ + β * ρₒᵤₜ
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
