"""
```julia
function lda(Z[, x]; conf, xc_func! = xc_lda(LDA_X, LDA_C_VWN), Vex = r -> -Z/r, δn=1e-8, maxit=100, xmax=25, Α=1e5)
```
Solve Kohn-Sham DFT Self-Consistent Field equations for an atom using
[local density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

On input:

* `Z`        - nuclear charge
* `x`        - grid, xᵢ = log(rᵢ), default is 401 points -35:0.1:5
* `conf`     - electronic configuration, i.e. `c"[He] 2s2 2p4"`, by default `conf = atom[Z]`
* `xc_func!` - LDA exchange-correlation functional, default is Slater exchange + VWN correlation
* `Vex`      - external potential, Coulomb -Z/r by default
* `δn`       - density convergence criterion
* `maxit`    - maximum permissible number of SCF iterations
* `xmax`     - grid is enlarged to this value for the Poisson equation solution
* `Α`        - parameter for the convex matrix pencil formulation

On output the function returns tuple: `(energy  = (total = E, kinetic = Ek, hartree = Eh, xc = Exc, potential = Ep), density = ρ, orbitals = ψ, iterations = Int(it))`
where

*  `lda(...).energy.total` - total energy, a.u.
*  `lda(...).density`      - electron density ρ(x) = 1/4π * ∑nᵢ*ψᵢ²(x)
*  `lda(...).orbitals`     - orbitals ψᵢ(x), energy levels εᵢ, azimuthal and radial quantum numbers l, nᵣ, and populations nᵢ (as listed in `conf`).

# Example
[Hooke atom](https://en.wikipedia.org/wiki/Hooke's_atom)

```jldoctest
julia> lda(2, Vex = r -> 1/8 * r^2).energy.total
[ Info: A grid of 401 points x = -35.0000:0.1000:5.0000 is used
[ Info: Starting SCF procedure with convergence threshold δn ≤ 1.000000e-08
[ Info: cycle           energy          δn
[ Info:   0           2.103849      1.89556590
[ Info:   1           2.013080      0.35778036
[ Info:   2           2.023521      0.03964354
[ Info:   3           2.026081      0.00360748
[ Info:   4           2.026212      0.00040386
[ Info:   5           2.026228      0.00004527
[ Info:   6           2.026229      0.00000518
[ Info:   7           2.026229      0.00000060
[ Info:   8           2.026229      0.00000007
[ Info:   9           2.026229      0.00000001
[ Info: SCF has converged after 9 iterations
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022579
│       EXCHANGE-CORRELATION          -0.523773
│       ELECTRON-NUCLEAR               0.899965
│       TOTAL ENERGY                   2.026229
└       VIRIAL RATIO                  -2.229260
2.0262293420299144
```
"""
function lda(Z::Real,
             x::AbstractRange = -35:0.1:5;             
          conf::Array{Array{T,1},1} = atom[Z],
      xc_func!::Function = xc_lda(LDA_X, LDA_C_VWN),
           Vex::Function = r -> -Z / r,
            δn::Real = 1e-8,
         maxit::Integer = 100,
          xmax::Real = 25,
             Α::Real = 1e5) where T <: Real
  
  μ = 1    # electron mass, a.u.
  β = 0.8  # initial value of the density mixing parameter

  # radial grid
  r, r², sqr, n, dx = exp.(x), exp.(2x), exp.(x/2), length(x), step(x)
  @info @sprintf("A grid of %3i points x = %5.4f:%5.4f:%5.4f is used", n, x[1], dx, x[end])

  ρᵢₙ, ρₒᵤₜ, prev_ρᵢₙ, prev_ρₒᵤₜ, V, vp, vh, vxc, εxc, ε = [similar(r) for _ = 1:10]

  Q = Nₑ(conf)

  # Thomas-Fermi density for neutral atom with Q electrons
  @. ρᵢₙ = TF(r, Q) * r
  @. vp  = Vex(r)

  # Indices of the main diagonal
  D = 1:n+1:n*n 

  H = -1/2μ*laplacian(x); K = H[D]
  S = copy(H)
  
  # Large grid for the Poisson equation
  xl = first(x):step(x):max(xmax, last(x))
  nl = length(xl)
  rmin, rmax = exp(first(xl)), exp(last(xl))
  
  # Differential operator of the Poisson equation
  Δ = laplacian(xl) - 1/4*I 
    
  @info @sprintf("Starting SCF procedure with convergence threshold δn ≤ %e", δn)
  @info "cycle\t\tenergy\t\tδn"

  Eₜₒₜ = 0.0; Δρ = 0.0; it = 0
  for i = 0:maxit
    # Poisson equation solution is Hartree potential, vh
    v = Δ \ vcat(-4π .* ρᵢₙ .* sqr .* r, zeros(nl - n))
    
    # Apply boundary conditions at r → 0 and r → ∞
    c₁ =  (√rmin * v[1] - √rmax * v[end] + Q * (1 - rmin)) / (rmax - rmin)
    c₂ = -(rmin  * √rmin * v[1] - √rmax * rmin * v[end] + Q * rmin * (1 - rmax)) / (rmax - rmin)
    @. @views vh = v[1:n] + c₁ * sqr + c₂ / sqr
    
    # Change variable vh(x) → vh(r)
    vh ./= sqr

    # Exchange-correlation potential and energy density
    xc_func!(ρᵢₙ ./ r, vxc, εxc)

    # Kohn-Sham potential
    @. V = (vp + vh + vxc) * r²

    # Schrödinger equation
    ∑nᵢεᵢ = 0.0; ρₒᵤₜ .= 0.0

    for l=0:length(conf)-1
      @. H[D] = K    + V + 1/2μ * (l + 1/2)^2
      @. S[D] = H[D] + Α * r²
      θ, y = eigen(Symmetric(H), Symmetric(S))
      @. ε = Α * θ / (1 - θ)
      for (nᵣ, nᵢ) in enumerate(conf[l+1])
        @views y[:, nᵣ] ./= sqrt(∫(dx, y[:, nᵣ] .^ 2 .* r²))
        @views ρₒᵤₜ .+= nᵢ / 4π * y[:, nᵣ] .^ 2
        ∑nᵢεᵢ += nᵢ * ε[nᵣ]
      end
    end

    # DFT total energy
    Eₜₒₜ = ∑nᵢεᵢ + 4π * ∫(dx, ρₒᵤₜ .* (-1/2 * vh .- vxc .+ εxc) .* r²)

    Δρ = 4π * ∫(dx, abs.(ρₒᵤₜ .- ρᵢₙ) .* r²)
    @info @sprintf "%3i\t%14.6f\t%14.8f" i Eₜₒₜ Δρ
    if Δρ < δn
      it = i
      break
    end

    # Convergence acceleration scheme
    if i > 0
      dρ = ∫(dx, ρₒᵤₜ .- prev_ρₒᵤₜ) / ∫(dx, ρᵢₙ .- prev_ρᵢₙ)
      β = min(max(1 / (1 - dρ), 0.1), 0.9)
    end
    @. prev_ρᵢₙ  = ρᵢₙ
    @. prev_ρₒᵤₜ = ρₒᵤₜ
    @. ρᵢₙ = (1 - β) * ρᵢₙ + β * ρₒᵤₜ
  end

  if Δρ > δn
    @warn "SCF does not converged after $maxit iterations!"
  else
    @info "SCF has converged after $it iterations"
  end

  # Wavefunctions (orbitals)
  @. V = vp + vh + vxc
  ∑ε, ρ, ψ = radial_shr_eq(V, x, conf = conf, μ = 1, Α = 1e5)
  
  # Total energy components:
  Ek  = ∑ε - 4π * ∫(dx, ρ .* V .* r²)
  Eh  = 2π * ∫(dx, ρ .* vh .* r²)
  Exc = 4π * ∫(dx, ρ .* εxc .* r²)
  Ep  = 4π * ∫(dx, ρ .* vp .* r²)
  E   = Ek + Ep + Eh + Exc
  
  # Virial ratio
  VR  = (Eh + Exc + Ep) / -Ek

  @info @sprintf("RESULTS SUMMARY:
      ELECTRON KINETIC       %16.6f
      ELECTRON-ELECTRON      %16.6f
      EXCHANGE-CORRELATION   %16.6f
      ELECTRON-NUCLEAR       %16.6f
      TOTAL ENERGY           %16.6f
      VIRIAL RATIO           %16.6f\n", Ek, Eh, Exc, Ep, E, VR)
      
  return (energy  = (total = E, kinetic = Ek, hartree = Eh, xc = Exc, potential = Ep), 
          density = ρ, orbitals = ψ, iterations = Int(it))
end
