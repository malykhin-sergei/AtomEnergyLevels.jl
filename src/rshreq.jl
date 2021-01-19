"""
    radial_shr_eq(V::Function = r -> -1/r, x::AbstractRange{T} = -30.0:0.1:5.0;
    conf = 1, μ = 1.0, α = 1e5) where {T}
    
    radial_shr_eq(V::Array{S}, x::AbstractRange{T}; 
    conf = 1, μ = 1.0, α = 1e5) where {S, T}

Solve Schrödinger equation ([atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are assumed)
```math
-\\frac{1}{2\\mu}\\left[\\frac{1}{r^{2}}\\frac{\\partial}{\\partial r}\\left(r^{2}\\frac{\\partial R_{nl}}{\\partial r}\\right)-\\frac{l(l+1)}{r^{2}}R_{nl}(r)\\right]+V(r)R_{nl}(r)=E_{nl}R_{nl}(r)
```
for [the particle in a spherically symmetric potential.](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)

## On input:

* `V` - potential, provided as an explicit function ``V(r)`` (by default `r -> -1/r`), or an array of values on a radial grid.
* `x` - uniform grid (`AbstractRange` object), such as ``r_i = \\exp(x_i)``. Default grid is `-30.0:0.1:5.0`
* `conf` - electronic configuration of interest, by default: ground state `conf = ((1.0),)`
* `μ` - [reduced mass](https://en.wikipedia.org/wiki/Reduced_mass), by default `μ = 1.0` (electron mass in a.u.)
* `α` - parameter ``α = 10^5`` is used in the matrix pencil ``\\mathrm{H} + α\\mathrm{S}`` for the generalized eigenproblem.

## On output function returns:

* sum of one-particle energies: ``E = \\sum_{i=1}^{N} n_i \\epsilon_i``
* particles density: ``\\rho(x) = \\frac{1}{4\\pi} \\sum_{i=1}^{N} n_i y_{i}^{2}(x)``
* orbitals ``y_i(x)``, corresponding eigenvalues ``\\epsilon_i`` (energy levels), azimuthal and radial quantum numbers ``l``, ``n_r``, and level populations ``n_i`` (as listed in `conf`).

## Examples
[Hydrogen atom](https://en.wikipedia.org/wiki/Hydrogen_atom)

```jldoctest
julia> radial_shr_eq(r -> -1/r, conf = conf_enc("1s1")).energy ≈ -1/2
true

julia> radial_shr_eq(r -> -1/r, conf = conf_enc("2s1")).energy ≈ -1/8
true

julia> radial_shr_eq(r -> -1/r, conf = conf_enc("2p1")).energy ≈ -1/8
true

julia> radial_shr_eq(r -> -1/r, conf = conf_enc("3s1")).energy ≈ -1/18
true

julia> radial_shr_eq(r -> -1/r, conf = conf_enc("3p1")).energy ≈ -1/18
true

julia> radial_shr_eq(r -> -1/r, conf = conf_enc("3d1")).energy ≈ -1/18
true
```
"""
function radial_shr_eq(V::AbstractArray, 
                       x::AbstractRange = -30.0:0.1:5.0; 
                    conf::Array{Array{T,1},1} = [[1]], 
                       μ::Real = 1, 
                       Α::Real = 1e5) where T <: Real

  if length(V) != length(x)
    throw(DimensionMismatch("V and r have different length"))
  end
  
  r, r², n, dx = exp.(x), exp.(2x), length(x), step(x)
  
  D = 1:n+1:n*n

  H = -1/2μ*laplacian(x); HD = H[D] .+ V .* r²
  S = copy(H)
 
  ε = similar(r); ρ = zero(r); ∑nᵢεᵢ = 0.0
  
  ψ = Dict()
  # the equation is solved separately for each subshell: s, p, d, f
  for (l, subshell) in enumerate(conf)
    @. H[D] = HD + 1/2μ * (l - 1/2)^2
    @. @views S[D] = H[D] + Α * r²
    θ, y = eigen(Symmetric(H), Symmetric(S))
    @. ε = Α * θ / (1 - θ)
    for (nᵣ, nᵢ) in enumerate(subshell)
      @views y[:, nᵣ] ./= sqrt(∫(dx, y[:, nᵣ] .^ 2 .* r²))
      @views ρ .+= nᵢ / 4π * y[:, nᵣ] .^ 2
      ∑nᵢεᵢ += nᵢ * ε[nᵣ]
      ψ[(nᵣ = nᵣ - 1, l = l - 1)] = (ϵᵢ = ε[nᵣ], nᵢ = nᵢ, ψᵢ = y[:, nᵣ])
    end
  end
  return (energy = ∑nᵢεᵢ, density = ρ, orbitals = ψ)
end
radial_shr_eq(V::Function = r -> -1/r, x::AbstractRange = -30.0:0.1:5.0; conf = 1, μ = 1, Α = 1e5) = 
              radial_shr_eq(V.(exp.(x)), x, conf = conf, μ = μ, Α = Α)
