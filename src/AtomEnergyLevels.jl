module AtomEnergyLevels

import LinearAlgebra: eigen!, diagind, Diagonal
import Printf: @sprintf

export laplacian, radial_shr_eq, TF, lda

include("dft_xc_functionals.jl")
export Slater, VWN

include("periodic_table.jl")
export atomic_shell, atomic_electron_configuration

"""
    laplacian(n, h)

Returns n×n differentiation matrix Δ = ∂²/∂x².

On input: `n` - size of a uniform grid with step size `h`.

# Example

Using the following code, one can find the vibrational levels of the radical OH⋅,
for which the potential energy surface is approximated through the Morse potential.

```https://en.wikipedia.org/wiki/Morse_potential```

Numerical solution is compared with analytic formula.
```julia-repl
julia> using LinearAlgebra, AtomEnergyLevels

julia> n = 100;

julia> D = 0.1994; β = 1.189; x₀ = 1.821;

julia> mH = 1.00794; mO = 15.9994; μ = 1822.8885*(mH*mO)/(mH+mO);

julia> xmin = 2.0 - x₀; xmax = 12.0 + x₀;

julia> x = range(xmin, xmax, length = n); dx = x.step.hi;

julia> V(x) = D*(exp(-β*(x-x₀))-1)^2 - D;  # Morse potential

julia> Δ = laplacian(n, dx);

julia> H = -1/2μ*Δ + Diagonal(V.(x));      # Hamiltonian

julia> ε, ψ = eigen(H);                    # numerical solution

julia> ω = β*sqrt(2D/μ); δ = ω^2 / 4D;

julia> E(n) = ω*(n+1/2) - δ*(n+1/2)^2 - D; # exact solution

julia> ε[1]
-0.19047200269535575

julia> E(0)
-0.19047201661032412

```
"""
function laplacian(n, h)
  Δ = Matrix{Float64}(undef, n, n)
  Δ[diagind(Δ, 0)] .= -1/3*π^2 / h^2
  for i=2:n
    Δ[diagind(Δ, i-1)] = Δ[diagind(Δ, 1-i)] .= 2*(-1)^i / (i-1)^2 / h^2
  end
  return Δ
end

∫(h, f) = h*sum(f) # rectangle rule for integration is accurate enough

"""
    radial_shr_eq(grid, V, conf, μ = 1.0, α = 1e5)

Solve Schrödinger equation (atomic units are assumed)

H = -1/2μ*Δ + V(r)

for [the particle in a spherically symmetric potential.](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)

On input:

* `grid` - logarithmic radial grid (a tuple (r, n, dx))
* `V`    - potential calculated on the grid points `r`
* `conf` - energy levels of interest.

For example, to find the energy levels with quantum numbers nᵣ = 1, 2, 3 and l = 0, 1, 2,
the levels of interest should be included in the tuple `conf` of the electronic
configuration as follows:

```
#          l=0        l=1        l=2
conf = ((1, 0, 0), (0, 0, 0), (0, 0, 0));
#   nᵣ = 0, 1, 2;   0, 1, 2;   0, 1, 2;
```

where the number nᵢ in a tuple is population of level with particular quantum numbers nᵣ, l.

As optional parameters, the particle mass `μ` and the coefficient `α` can be specified.
Parameter `α` is used in the matrix pencil `H + α*Diagonal(r²)` for the generalized
eigenproblem solution.

On output:

* sum of the one-particle energies: E = ∑nᵢ*εᵢ
* particles density: ρ = 1/4π * ∑nᵢ*ψᵢ²(r)
* matrix, containing complete set of orbitals (values of the wavefunctions ψᵢ(r) on the
  radial grid `r`, corresponding eigenvalues εᵢ (energy levels), azimuthal and
  radial quantum numbers l, nᵣ, and level populations nᵢ (as listed in `conf`).

# Examples
## 3D isotropic harmonic oscillator problem

```https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator```

Exact eigenvalues are: E = ħω*(2nᵣ + l + 3/2).
The numerical solution is following.

```julia-repl
julia> r, n, dx = begin
           x = range(-30, 20, length = 501)
           exp.(x), x.len, x.step.hi
       end;

julia> conf = ((1,0,0,), (0,0,0), (0,0,0));

julia> ψ = radial_shr_eq((r, n, dx), 1/2*r.^2, conf).orbitals;

julia> ψ[1:3,:]
3×9 Array{Float64,2}:
 0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0  2.0  # azimuthal quantum numbers, l
 0.0  1.0  2.0  0.0  1.0  2.0  0.0  1.0  2.0  # radial quantum numbers, nᵣ
 1.5  3.5  5.5  2.5  4.5  6.5  3.5  5.5  7.5  # energy levels
```

## Hydrogen atom

```https://en.wikipedia.org/wiki/Hydrogen_atom```

Hamiltonian of the particle is:

H = -1/2*Δ - 1/r

and exact solution is:

E = -1/2(nᵣ + l + 1)², R1s(r) = 2*exp(-r)

```julia-repl
julia> r, n, dx = begin
           x = range(-30, 20, length = 501)
           exp.(x), x.len, x.step.hi
       end;

julia> result = radial_shr_eq((r, n, dx), -1 ./ r, 1);

julia> result.energy
-0.5000000000310655

julia> ρ = result.density;

julia> 4π * dx * sum(ρ .* r .^3)
0.9999999999999998

julia> R1s = -result.orbitals[5:end,1];

julia> R1s_exact = 2*exp.(-r);

julia> dx * sum((R1s .- R1s_exact).^2 .* r .^3)
9.313469540251104e-22
```
"""
function radial_shr_eq((r, n, dx), V, conf; μ = 1.0, α = 1e5)
  r² = r .* r

  H = -1/2μ*laplacian(n, dx) + Diagonal(V .* r²)

  ρ = zeros(n)
  ψ = zeros(n + 4, length(collect(Iterators.flatten(conf))))

  k = 1; E = 0.0

  for (l, subshell) in enumerate(conf)

    Hl = H + Diagonal(fill(1/2μ * (l - 1/2)^2, n))
    θ, y = eigen!(Hl, Hl + Diagonal(α * r²))
    ε = α*θ ./ (1 .- θ)

    for (i, occ) in enumerate(subshell)
      # normalize wavefunctions
      y[:, i] *= 1.0 / sqrt(∫(dx, y[:, i] .^ 2 .* r²))

      # make density from occupied orbitals
      ρ .+= occ / 4π * y[:, i] .^ 2 ./ r

      # calculate total energy
      E += occ * ε[i]

      # save results
      ψ[1, k] = l-1                       # azimuthal quantum number
      ψ[2, k] = i-1                       # radial quantum number
      ψ[3, k] = ε[i]                      # level energy
      ψ[4, k] = occ                       # level occupation
      ψ[5:end, k] = y[:, i] ./ sqrt.(r)   # orbital
      k += 1
    end
  end
  return (energy = E, density = ρ, orbitals = ψ)
end

"""
    TF(r, Z)

Approximate solution to the Thomas-Fermi equation for neutral atom with nuclear charge `Z`.
Returns density at the distance `r`.

```https://en.wikipedia.org/wiki/Thomas-Fermi_model```

Moliere, G. (1947). Theorie der streuung schneller geladener teilchen
i. einzelstreuung am abgeschirmten coulomb-feld.
Zeitschrift für Naturforschung A, 2(3), 133-145.

# Example
```julia-repl
julia> r, n, dx = begin
           x = range(-30, 20, length = 501)
           exp.(x), x.len, x.step.hi
       end;

julia> 4π * dx * sum(TF.(r, 86) .* r .^3)
85.99999999999511
```
"""
function TF(r, Z)
  b = 1/2*(3π/4)^(2/3)*Z^(-1/3); B = [0.1 0.55 0.35]; β = [6.0 1.20 0.30]
  return Z/4π/b^2/r*(B[1]*β[1]^2*exp(-β[1]*r/b) +
                     B[2]*β[2]^2*exp(-β[2]*r/b) +
                     B[3]*β[3]^2*exp(-β[3]*r/b))
end

"""
```julia
function lda((r, n, dx) = begin
                x = range(-30, 20, length = 501);
                    exp.(x), x.len, x.step.hi end;
                Z = nothing,
             conf = nothing,
               xc = SVWN,
               vp = nothing,
             ρ_in = nothing,
                β = 0.3,
               δn = 1.0e-6,
               δE = 5.0e-7,
            maxit = 100,
                μ = 1,
                α = 1e5)
```
Solve Kohn-Sham DFT Self-Consistent Field equations for an atom using
[local density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

On input:

* `(r, n, dx)` - radial grid, dr = r * dx
   (default is r::501-element Array{Float64,1}, n = 501, dx = 0.1)
* `Z`      - nuclear charge (optional)
* `conf`   - electronic configuration
   (optional, but `Z` and/or `conf` must be provided anyway)
* `xc`     - LDA exchange-correlation functional
   (default is S.H. Vosko, L. Wilk, and M. Nusair functional, SVWN(ρ))
*  `vp`    - external potential
   (default is Coulomb potential `vp = -Z/r`)
*  `ρ_in`  - input density
   (default is Thomas-Fermi density)
*  `β`     - density mixing parameter
   (default is `β = 0.3` which almost guarantees convergence)
*  `δn`    - density convergence criterion
   (default is `δn = 1.0e-6`)
*  `δE`    - energy convergence criterion
   (default is `δE = 5.0e-7`)
*  `maxit` - maximum permissible number of SCF iterations
   (optional, default is `maxit = 100`, enough for most tasks)

On output:

*  `lda(...).energy`   - total energy, a.u.
*  `lda(...).density`  - electron density ρ(r) = 1/4π * ∑nᵢ*ψᵢ²(r)
*  `lda(...).orbitals` - matrix, containing complete set of orbitals
    (values of the wavefunctions ψᵢ(r) on the radial grid `r`,
    corresponding eigenvalues εᵢ (energy levels), azimuthal and
    radial quantum numbers l, nᵣ, and level populations nᵢ (as listed in `conf`).

# Example
## Helium
```julia-repl
julia> lda(Z = 2, β = 0.8).energy;
[ Info: Neutral atom with Z =  2 and He electron configuration is assumed.
[ Info: Using logarithmic 501 point grid and step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -3.231952	    1.517742
[ Info:   1	     -2.838866	    0.275522
[ Info:   2	     -2.843659	    0.048204
[ Info:   3	     -2.835333	    0.009764
[ Info:   4	     -2.835105	    0.001852
[ Info:   5	     -2.834863	    0.000380
[ Info:   6	     -2.834844	    0.000076
[ Info:   7	     -2.834837	    0.000016
[ Info:   8	     -2.834836	    0.000003
[ Info:   9	     -2.834836	    0.000001
[ Info:  10	     -2.834836	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               2.767922
│       ELECTRON-ELECTRON              1.996120
│       EXCHANGE-CORRELATION          -0.973314
│       ELECTRON-NUCLEAR              -6.625564
│       TOTAL ENERGY                  -2.834836
└       VIRIAL RATIO                   2.024175
```
## Lithium atom EA and IP
Electron affinity and ionization potential for an atom can be calculated
using Slater-Janak transition state.
```julia-repl
julia> const ev = 27.21138400
27.211384

julia> EA = -lda(Z = 3, conf = ((2, 1.5),()), β = 0.8).orbitals[3,end]*ev
[ Info: Using logarithmic 501 point grid and step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -6.539026	    3.619811
[ Info:   1	     -7.959516	    2.552724
[ Info:   2	     -7.352136	    0.792187
[ Info:   3	     -7.361211	    0.205290
[ Info:   4	     -7.367748	    0.036280
[ Info:   5	     -7.367463	    0.007823
[ Info:   6	     -7.367563	    0.001581
[ Info:   7	     -7.367568	    0.000355
[ Info:   8	     -7.367572	    0.000081
[ Info:   9	     -7.367573	    0.000020
[ Info:  10	     -7.367573	    0.000005
[ Info:  11	     -7.367573	    0.000001
[ Info:  12	     -7.367573	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               7.264385
│       ELECTRON-ELECTRON              4.323666
│       EXCHANGE-CORRELATION          -1.707320
│       ELECTRON-NUCLEAR             -17.248303
│       TOTAL ENERGY                  -7.367573
└       VIRIAL RATIO                   2.014205
0.7269491150354416

julia> IP = -lda(Z = 3, conf = ((2, 0.5),()), β = 0.8).orbitals[3,end]*ev
[ Info: Using logarithmic 501 point grid and step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -8.080241	    1.621832
[ Info:   1	     -7.297013	    0.316934
[ Info:   2	     -7.261408	    0.082710
[ Info:   3	     -7.259677	    0.021497
[ Info:   4	     -7.259888	    0.005607
[ Info:   5	     -7.260061	    0.001463
[ Info:   6	     -7.260127	    0.000381
[ Info:   7	     -7.260148	    0.000099
[ Info:   8	     -7.260154	    0.000026
[ Info:   9	     -7.260155	    0.000007
[ Info:  10	     -7.260156	    0.000002
[ Info:  11	     -7.260156	    0.000000
[ Info:  12	     -7.260156	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               7.168653
│       ELECTRON-ELECTRON              3.619145
│       EXCHANGE-CORRELATION          -1.588566
│       ELECTRON-NUCLEAR             -16.459388
│       TOTAL ENERGY                  -7.260156
└       VIRIAL RATIO                   2.012764
5.316382601516975
```
Experimental values are: IP = 5.39 eV, EA = 0.62 eV.

## Hooke atom
`https://en.wikipedia.org/wiki/Hooke's_atom`

An exact solution for Hooke's atom is E = 2.0 a.u.

```julia-repl
julia> r, n, dx = begin
    x = range(-30, 20, length = 501)
    exp.(x), x.len, x.step.hi
end;

julia> lda(conf = 2, vp = 1/8 * r.^2, β = 0.8).energy;
[ Info: Neutral atom with Z =  2 is assumed.
[ Info: Using logarithmic 501 point grid and step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	      2.103849	    1.895566
[ Info:   1	      2.013080	    0.357780
[ Info:   2	      2.022087	    0.073897
[ Info:   3	      2.025397	    0.014779
[ Info:   4	      2.026076	    0.003000
[ Info:   5	      2.026201	    0.000621
[ Info:   6	      2.026224	    0.000130
[ Info:   7	      2.026228	    0.000028
[ Info:   8	      2.026229	    0.000006
[ Info:   9	      2.026229	    0.000001
[ Info:  10	      2.026229	    0.000000
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
function lda((r, n, dx) = begin
                x = range(-30, 20, length = 501);
                    exp.(x), x.len, x.step.hi end;
                Z = nothing,
             conf = nothing,
               xc = ρ -> Slater(ρ) .+ VWN(ρ),
               vp = nothing,
             ρ_in = nothing,
                β = 0.3,
               δn = 1.0e-6,
               δE = 5.0e-7,
            maxit = 100,
                μ = 1,
                α = 1e5)

  # Check input parameters
  if Z == nothing && conf == nothing
    @error "System is not specified: provide ether atomic number Z or/and electron configuration."
    return (energy = nothing, density = nothing, orbitals = nothing)
  elseif Z == nothing
    Q = sum(Iterators.flatten(conf))
    Z = Q
    @info @sprintf("Neutral atom with Z = %2i is assumed.", Z)
  elseif conf == nothing
    conf = atomic_electron_configuration[Z]
    Q = Z
    s = keys(atomic_electron_configuration)[Z]
    @info @sprintf("Neutral atom with Z = %2i and %s electron configuration is assumed.", Z, s)
  else
    Q = sum(Iterators.flatten(conf))
  end

  r² = r .* r; r³ = r² .* r; sqr = sqrt.(r)
  @info @sprintf("Using logarithmic %3i point grid with step dx = %5.4f", n, dx)

  if vp == nothing vp = -Z ./ r end

  if ρ_in == nothing
    ρ_in = TF.(r, Q)
    @info "Using Thomas-Fermi starting electron density"
  end

  nlevels = length(collect(Iterators.flatten(conf)))
  ψ = zeros(n + 4, nlevels)

  Δ = laplacian(n, dx)
  L = Δ - Diagonal(fill(1/4, n))

  Etot = 0.0; Δρ = 0.0
  ρ_out, V, vh, vxc, εxc = [Array{Float64,1}(undef, n) for _ = 1:5]

  @info @sprintf("Starting SCF procedure with density mixing parameter β = %5.4f
      and convergence threshold |Δρ| ≤ %e", β, δn)
  @info "cycle\t\tenergy\t\t|Δρ|"
  for i = 0:maxit
    # Solve the Poisson equation to obtain Hartree potential
    vh = L \ (-4π * ρ_in .* sqr .* r²)
    # Apply boundary conditions at r → 0 and r → ∞
    vh .-= (vh[n] - Q / sqr[n])/sqr[n] .* sqr .+
           (vh[1] - Q * sqr[1])*sqr[1] ./ sqr
    # Change variable vh(x) → vh(r)
    vh ./= sqr

    # Calculate exchange-correlation potential and energy density
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#the-xc-term)
    for i=1:n vxc[i], εxc[i] = xc(ρ_in[i]) end

    # Assemble Kohn-Sham potential
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#kohn-sham-equations)
    V = vp + vh + vxc

    # Solve the Schrödinger equation to find new density - ρ_out
    # and bands energy ∑ε
    H = -1/2μ*Δ + Diagonal(V .* r²)

    ∑ε = 0.0; ρ_out .= 0.0

    # the equation is solved separately for each subshell: s, p, d, f
    for (l, subshell) in enumerate(conf)
      Hl = H + Diagonal(fill(1/2μ * (l - 1/2)^2, n))
      θ, y = eigen!(Hl, Hl + α*Diagonal(r²))
      ε = α*θ ./ (1 .- θ)
      for (i, occ) in enumerate(subshell)
        # normalize wavefunctions
        y[:, i] *= 1 / sqrt(∫(dx, y[:, i] .^ 2 .* r²))
        # build density
        ρ_out .+= occ / 4π * y[:, i] .^ 2 ./ r
        # sum up the energy levels
        ∑ε += occ * ε[i]
      end
    end

    E_prev = Etot

    # DFT total energy:
    # (https://www.theoretical-physics.net/dev/quantum/dft.html#total-energy)
    Etot = ∑ε + 4π * ∫(dx, ρ_out .* (-1/2 * vh .- vxc .+ εxc) .* r³)

    Δρ = 4π * ∫(dx, abs.(ρ_out - ρ_in) .* r³)
    @info @sprintf "%3i\t%14.6f\t%12.6f\n" i Etot Δρ

    # density converged if there are no more charge oscillations
    Δρ < δn && abs(Etot - E_prev) < δE && break

    # admix density from previous iteration to converge SCF
    ρ_in = (1.0 - β) * ρ_in + β * ρ_out
  end
  Δρ > δn && @warn "SCF did not converge after $maxit iterations!"

  # obtain orbitals
  ∑ε, ρ, ψ = radial_shr_eq((r, n, dx), V, conf)

  # remove numerical noise
  ρ = [ifelse(x > eps(), x, 0.0) for x in ρ]

  # total energy components
  Ek  = ∑ε - 4π * ∫(dx, ρ .* V .* r³)
  Eh  = 2π * ∫(dx, ρ .* vh .* r³)
  Exc = 4π * ∫(dx, ρ .* εxc .* r³)
  Ep  = 4π * ∫(dx, ρ .* vp .* r³)
  E   = Ek + Ep + Eh + Exc
  VR  = (Eh + Exc + Ep) / -Ek

  @info @sprintf("RESULTS SUMMARY:
      ELECTRON KINETIC       %16.6f
      ELECTRON-ELECTRON      %16.6f
      EXCHANGE-CORRELATION   %16.6f
      ELECTRON-NUCLEAR       %16.6f
      TOTAL ENERGY           %16.6f
      VIRIAL RATIO           %16.6f\n", Ek, Eh, Exc, Ep, E, VR)
  return (energy = E, density = ρ, orbitals = ψ)
end
end # module
