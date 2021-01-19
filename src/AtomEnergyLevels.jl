module AtomEnergyLevels

include("integrate.jl")
include("sincdif.jl")
include("conf_parse.jl")
include("periodic_table.jl")
include("dft_xc_functionals.jl")
include("rshreq.jl")
include("TF.jl")

import LinearAlgebra: I, eigen, diagind, Symmetric
import Printf: @sprintf

export laplacian, radial_shr_eq, TF, lda
export LDA_X, LDA_C_CHACHIYO, LDA_C_VWN
export SVWN!, Xα!
export conf_enc, Nₑ, @c_str
export shells, atom
export laplacian
export radial_shr_eq

"""
```julia
function lda(Z, x = -30.0:0.1:20.0; conf = atom[Z],
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

```@example
lda(2, Vex = r -> 1/8 * r^2).energy.total;
```
"""
function lda(Z::Real,
             x::AbstractRange = -35:0.1:5;             
          conf = atom[Z],
           xc!::Function = SVWN!,
           Vex::Function = r -> -Z / r,
            δn::Real = 1e-8,
         maxit::Integer = 100,
             μ::Real = 1,
          xmax::Real = 25,
             Α::Real = 1e5)

  β = 0.8  # initial value of the density mixing parameter

  # radial grid
  r, r², sqr, n, dx = exp.(x), exp.(2x), exp.(x/2), length(x), step(x)
  @info @sprintf("A grid of %3i points x = %5.4f:%5.4f:%5.4f is used", n, x[1], dx, x[end])

  ρᵢₙ, ρₒᵤₜ, prev_ρᵢₙ, prev_ρₒᵤₜ, V, vp, vh, vxc, εxc, ε = [similar(r) for _ = 1:10]

  Q = Nₑ(conf)
  @. ρᵢₙ = TF(r, Q) * r
  @. vp  = Vex(r)

  D = 1:n+1:n*n

  H = -1/2μ*laplacian(x); HD = H[D]
  S = copy(H)
  
  xl = first(x):step(x):xmax
  nl = length(xl)
  rmin, rmax = exp(first(xl)), exp(last(xl)) 
  Δ = laplacian(xl) - 1/4*I 
    
  Eₜₒₜ = 0.0; Δρ = 0.0;  

  @info @sprintf("Starting SCF procedure with convergence threshold |Δρ| ≤ %e", δn)
  @info "cycle\t\tenergy\t\t|Δρ|"
  for i = 0:maxit
    # Solve the Poisson equation to obtain Hartree potential
    v = Δ \ vcat(-4π .* ρᵢₙ .* sqr .* r, zeros(nl - n))
    # Apply boundary conditions at r → 0 and r → ∞
    c₁ =  (√rmin * v[1] - √rmax * v[end] + Q * (1 - rmin)) / (rmax - rmin)
    c₂ = -(rmin  * √rmin * v[1] - √rmax * rmin * v[end] + Q * rmin * (1 - rmax)) / (rmax - rmin)
    @. vh = @views v[1:n] + c₁ * sqr + c₂ / sqr
    # Change variable vh(x) → vh(r)
    @. vh /= sqr

    # Exchange-correlation potential and energy density
    xc!(ρᵢₙ ./ r, vxc, εxc)

    # Kohn-Sham potential
    @. V = (vp + vh + vxc) * r²

    # Solve the Schrödinger equation to find new density
    # and bands energy ∑nᵢεᵢ

    ∑nᵢεᵢ = 0.0
    ρₒᵤₜ .= 0.0

    # the equation is solved separately for each subshell: s, p, d, f
    for (l, subshell) in enumerate(conf)
      @. H[D] = HD   + V + 1/2μ * (l - 1/2)^2
      @. S[D] = H[D] + Α * r²
      θ, y = eigen(Symmetric(H), Symmetric(S))
      @. ε = Α * θ / (1 - θ)
      for (nᵣ, nᵢ) in enumerate(subshell)
        @views y[:, nᵣ] /= sqrt(∫(dx, y[:, nᵣ] .^ 2 .* r²))
        @views ρₒᵤₜ .+= nᵢ / 4π * y[:, nᵣ] .^ 2
        ∑nᵢεᵢ += nᵢ * ε[nᵣ]
      end
    end
    E = Eₜₒₜ
    # DFT total energy:
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#total-energy)
    Eₜₒₜ = ∑nᵢεᵢ + 4π * ∫(dx, ρₒᵤₜ .* (-1/2 * vh .- vxc .+ εxc) .* r²)

    Δρ = 4π * ∫(dx, abs.(ρₒᵤₜ .- ρᵢₙ) .* r²)
    @info @sprintf "%3i\t%14.6f\t%12.6f" i Eₜₒₜ Δρ

    # density converged if there are no more charge oscillations
    Δρ < δn && break

    if i > 0
      dρₒᵤₜ = ∫(dx, (ρₒᵤₜ .- prev_ρₒᵤₜ) .* r) / ∫(dx, (ρᵢₙ .- prev_ρᵢₙ) .* r)
      β = min(max(1/(1 - dρₒᵤₜ), 0.1), 0.9)
    end
    @. prev_ρᵢₙ  = ρᵢₙ
    @. prev_ρₒᵤₜ = ρₒᵤₜ
    # admix density from previous iteration to converge SCF
    @. ρᵢₙ = (1 - β) * ρᵢₙ + β * ρₒᵤₜ
  end
  Δρ > δn && @warn "SCF does not converged after $maxit iterations!"

  # get orbitals
  @. V = vp + vh + vxc
  ∑ε, ρ, ψ = radial_shr_eq(V, x, conf = conf, μ = μ, Α = 1e5)
  
  # total energy components:
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
