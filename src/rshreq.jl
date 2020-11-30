"""
    radial_shr_eq(V, x, conf, μ, α)

Solve Schrödinger equation (atomic units are assumed)

    H = -1/2μ*Δ + V(r)

for [the particle in a spherically symmetric potential.](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)

# On input:

* `V` - potential, provided as explicit function V(r) (by default `r -> -1/r`), or array of values on a radial grid r.
* `x` - uniform grid (range object), such as rᵢ = exp(xᵢ). Default grid is `-30.0:0.1:5.0`
* `conf` - electronic configuration of interest, by default: ground state `conf = ((1.0),)`
* `μ` - reduced mass, by default `μ = 1.0` (electron mass in a.u.)
* `α` - parameter `α = 1e5` is used in the matrix pencil `H + α*Diagonal(r²)` for the generalized eigenproblem solution.

# On output function returns:

* sum of one-particle energies: E = ∑nᵢ*εᵢ
* particles density: ρ(x) = 1/4π * ∑nᵢ*ψᵢ²(x)
* orbitals ψᵢ(x), corresponding eigenvalues εᵢ (energy levels), azimuthal and radial quantum numbers l, nᵣ, and level populations nᵢ (as listed in `conf`).

# Examples
[Hydrogen atom](https://en.wikipedia.org/wiki/Hydrogen_atom)

```julia-repl
julia> radial_shr_eq(r -> -1/r).energy # n = 1 (ground state)
-1.9999999999812186

julia> radial_shr_eq(r -> -1/r, conf = ((0, 1),)).energy # n = 2 (excited state)
-0.12499999998289188
```
"""
radial_shr_eq(V::Function = r -> -1/r, x::AbstractRange{T} = -30.0:0.1:5.0;
              conf = 1, μ = 1.0, α = 1e5) where {T} = 
              _radial_shr_eq(x, V.(exp.(x)), conf, μ, α)

radial_shr_eq(V::Array{S}, x::AbstractRange{T}; 
              conf = 1, μ = 1.0, α = 1e5) where {S, T} = 
              _radial_shr_eq(x, V, conf, μ, α)

function _radial_shr_eq(x, V, conf, μ, α)
  if length(V) != length(x)
    throw(DomainError(x, "Arrays V and r have different length"))
  end
  
  r, r², n, dx = exp.(x), exp.(2x), length(x), step(x)
  
  H = -1/2μ*laplacian(n, dx) + Diagonal(V .* r²)

  E = 0.0; ρ = zeros(n); ψ = Dict()
  
  for (l, subshell) in enumerate(conf)
    Hl = H + Diagonal(fill(1/2μ * (l - 1/2)^2, n))
    θ, y = eigen!(Hl, Hl + Diagonal(α * r²))
    ε = α*θ ./ (1 .- θ)
    for (nᵣ, occ) in enumerate(subshell)
      y[:, nᵣ] /= sqrt(∫(dx, y[:, nᵣ] .^ 2 .* r²))
      ρ .+= occ / 4π * y[:, nᵣ] .^ 2
      E += occ * ε[nᵣ]
      ψ[(nᵣ = nᵣ - 1, l = l - 1)] = (ϵᵢ = ε[nᵣ], nᵢ = occ, ψᵢ = y[:, nᵣ])
    end
  end
  return (energy = E, density = ρ, orbitals = ψ)
end
