module AtomEnergyLevels

import LinearAlgebra: eigen, diag, diagind
import Printf: @sprintf

export radial_grid, laplacian, radial_shr_eq, TF, lda

include("dft_xc_functionals.jl")
export Xα, SVWN

include("periodic_table.jl")
export atomic_shell, atomic_electron_configuration

"""
    radial_grid(n = 401, xmin = -30, xmax = 20)

Given a size `n`, computes a logarithmic grid on the radial coordinate
from `r = exp(xmin)` to `exp(xmax)`, returning an object which contains 
grid values (vector of size `n`), `n`, and step dx = dr / r.

# Example
```julia-repl
julia> r, n, dx = radial_grid(51);

julia> r[end]-r[1]*exp((n-1)*dx)
0.0

julia> n
51
```
"""
function radial_grid(n = 401, xmin = -30, xmax = 20)
  x = range(xmin, xmax, length = n)
  return (r = exp.(x), n = x.len, dx = x.step.hi)
end


"""
    laplacian(n, h)
    
Create matrix representation of the Laplace operator Δ = ∂²/∂x² 
for the 1D Schrödinger equation (atomic units are assumed).

H = -1/2μ*Δ + V(x)
    
On input: `n` - size of grid, `h` - grid step.
On output: n×n differentiation matrix Δ = ∂²/∂x².

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

# simple integration via rectangle rule
∫(h, f) = h*sum(f)

"""
    radial_shr_eq(grid, V, conf, μ = 1.0)
    
Solve Schrödinger equation (atomic units are assumed)

H = -1/2μ*Δ + V(r)

for [the particle in a spherically symmetric potential.](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)
Logarithmic radial `grid` created with `log_radial_grid()` function, potential `V` 
calculated on the grid points and energy levels of interest (`conf` variable) should be 
provided on input. As an optional parameter, the particle mass `μ` can be specified. 

The function returns: sum of the one-particle energies, 

E = ∑εᵢ

```julia
    E = radial_shr_eq(grid, V, conf, μ).energy
```

the probability density

ρ = ∑nᵢψ²ᵢ

```julia
    ρ = radial_shr_eq(grid, V, conf, μ).density
```
 
according to the population nᵢ of levels listed in the `conf` tuple, 
and a table of orbitals (values of the wavefunction ψᵢ(r) on the radial grid) 
with corresponding eigenvalues εᵢ (energy levels), azimuthal and 
radial quantum numbers (l, nᵣ), and populations of energy levels nᵢ.

```julia
    ψ = radial_shr_eq(grid, V, conf, μ).orbitals
```

!!! note 

    `log_radial_grid(...).orbitals` actually returns wavefunctions yᵢ(x), 
    where x - uniform grid. The radial grid is `r = exp.(x)`. 
    True radial wavefunction can be obtained by transformation `R = y ./ sqrt.(r)`.

# Examples
## 3D isotropic harmonic oscillator problem

```https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator```

To find the energy levels with quantum numbers nᵣ = 1, 2, 3 and l = 0, 1, 2, 
the levels of interest should be included in the tuple of the electronic configuration 
as follows

```julia
#          l=0        l=1        l=2
conf = ((0, 0, 0), (0, 0, 0), (0, 0, 0));
#   nᵣ = 0, 1, 2;   0, 1, 2;   0, 1, 2;  
```
where the tuples `(0,0,0)` indicate population numbers of each shell.

The output should be a matrix

```julia
3×9 Array{Float64,2}:
 0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0  2.0  # azimuthal quantum numbers, l
 0.0  1.0  2.0  0.0  1.0  2.0  0.0  1.0  2.0  # radial quantum numbers, nᵣ
 1.5  3.5  5.5  2.5  4.5  6.5  3.5  5.5  7.5  # energy levels
```
Exact eigenvalues are: E = ħω*(2nᵣ + l + 3/2). The numerical solution
is:
```julia-repl
julia> r, n, dx = radial_grid(301, -20, 5);

julia> conf = ((0,0,0,), (0,0,0), (0,0,0));

julia> ψ = radial_shr_eq((r, n, dx), 1/2*r.^2, conf).orbitals;

julia> ψ[1:3,:]
3×9 Array{Float64,2}:
 0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0  2.0
 0.0  1.0  2.0  0.0  1.0  2.0  0.0  1.0  2.0
 1.5  3.5  5.5  2.5  4.5  6.5  3.5  5.5  7.5
```
## Hydrogen-like atom

```https://en.wikipedia.org/wiki/Hydrogen-like_atom```

Hamiltonian of the particle is:  

H = -1/2μ*Δ - Z/r

and exact solution is:

E = -μ*Z²/2n²

```julia-repl
julia> (r, n, dx) = radial_grid(500, -30, 10);

julia> ψ = radial_shr_eq((r, n, dx), -2 ./ r, 1, 0.25).orbitals;

julia> ψ[1:3,:]
3×1 Array{Float64,2}:
  0.0
  0.0
 -0.5000000000007793
```
"""
function radial_shr_eq(grid, V, conf, μ = 1.0)
  r, n, dx = grid
  r² = r .* r
  
  H = -1/2μ*laplacian(n, dx)
  for i=1:n H[i,:] ./= r²[i] end
  T = diag(H);
  
  ρ = zeros(n)
  ψ = zeros(n + 4, length(collect(Iterators.flatten(conf))))
  
  k = 1; E = 0.0

  for (l, subshell) in enumerate(conf)

    H[diagind(H, 0)] = T + V .+ 1/2μ * (l - 1/2)^2 ./ r²
    ε, y = eigen(H, sortby=real)

    for (i, occ) in enumerate(subshell)
      # normalize wavefunctions
      y[:, i] *= 1.0 / sqrt(∫(dx, y[:, i] .^ 2 .* r²))

      # make density from occupied orbitals
      ρ .+= occ / 4π * y[:, i] .^ 2 ./ r
      
      # calculate total energy
      E += occ * ε[i]

      # write results to the matrix
      ψ[1, k] = l-1          # azimuthal quantum number
      ψ[2, k] = i-1          # radial quantum number
      ψ[3, k] = ε[i]         # energy level
      ψ[4, k] = occ          # occupation
      ψ[5:end, k] = y[:, i]  # orbital
      k += 1
    end
  end
  return (energy = E, density = ρ, orbitals = ψ)
end

"""
    `TF(r, Z)`
    
Approximate solution to the Thomas-Fermi equation for neutral atom with nuclear charge `Z`.
Returns density at the distance `r`.

```https://en.wikipedia.org/wiki/Thomas-Fermi_model```

Moliere, G. (1947). Theorie der streuung schneller geladener teilchen 
i. einzelstreuung am abgeschirmten coulomb-feld. 
Zeitschrift für Naturforschung A, 2(3), 133-145.

# Example
```julia-repl
julia> (r, n, dx) = radial_grid(400, -30, 20); 

julia> 4π * dx * sum(TF.(r, 86) .* r .^3)
85.99999999999687

```
"""
function TF(r, Z)
  b = 1/2*(3π/4)^(2/3)*Z^(-1/3); B = [0.1 0.55 0.35]; β = [6.0 1.20 0.30]
  return Z/4π/b^2/r*(B[1]*β[1]^2*exp(-β[1]*r/b) +
                     B[2]*β[2]^2*exp(-β[2]*r/b) +
                     B[3]*β[3]^2*exp(-β[3]*r/b))
end

function radial_shr_eq!(grid, V, conf, H, ρ, ψ, μ = 1.0)
  r, n, dx = grid
  r² = r .* r
  T = diag(H)
  k = 1; E = 0.0; ρ .= 0.0
  for (l, subshell) in enumerate(conf)
    H[diagind(H, 0)] = T + V .+ 1/2μ * (l - 1/2)^2 ./ r²
    ε, y = eigen(H, sortby=real)
    for (i, occ) in enumerate(subshell)
      y[:, i] *= 1 / sqrt(∫(dx, y[:, i] .^ 2 .* r²))
      ρ .+= occ / 4π * y[:, i] .^ 2 ./ r
      E += occ * ε[i]      
      ψ[1, k] = l-1
      ψ[2, k] = i-1
      ψ[3, k] = ε[i]
      ψ[4, k] = occ
      ψ[5:end, k] = y[:, i]
      k += 1
    end
  end
  return E
end

function denoise!(arr, ξ = eps())
  for i=1:length(arr)
    if arr[i] < ξ 
       arr[i] = 0.0 
    end
  end
end

"""
```julia
    lda(grid = radial_grid(401, -30, 20);
           Z = nothing, 
        conf = nothing,
          xc = SVWN,
          vp = nothing,
        ρ_in = nothing,
           β = 0.3,
          δn = 1.0e-6,
          δE = 0.5e-6
       maxit = 50)
```

Solve DFT Kohn-Sham Self-Consistent Field equations for an atom.

On input:

* `grid`   - radial grid
* `Z`      - nuclear charge (optional)
* `conf`   - electronic configuration 
   (optional, but `Z` and/or `conf` must be provided anyway)
* `xc`     - LDA exchange-correlation functional 
   (some function f(ρ), default is S.H. Vosko, L. Wilk, and M. Nusair functional, SVWN)
*  `vp`    - external potential 
   (default is Coulomb potential `vp = -Z/r`)
*  `ρ_in`  - initial density 
   (default is Thomas-Fermi density)
*  `β`     - density mixing parameter
   (default is `β = 0.3` which almost guarantees convergence)
*  `δn`    - density convergence criterion
   (default is `δn = 1.0e-6`)
*  `δE`    - energy convergence criterion 
   (default is `δE = 5.0e-7`)     
*  `maxit` - maximum permissible number of SCF iterations
   (optional, default is `maxit = 50`, enough for most tasks)

On output function returns:

*  `lda(...).energy`   - total energy, a.u.
*  `lda(...).density`  - density ρ(r), array size `n`, where
    `n` - radial grid size
*  `lda(...).orbitals` - table of orbitals (values of the 
    wavefunction ψᵢ(r) on the radial grid) with corresponding 
    eigenvalues εᵢ (energy levels), azimuthal and radial quantum 
    numbers (l, nᵣ), and populations nᵢ.

# Example
## Helium
```julia-repl
julia> lda(Z = 2, β = 0.8).energy;
[ Info: Neutral atom with Z =  2 and He electron configuration is assumed.
[ Info: Using logarithmic 401 point grid and step dx = 0.1250
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000 
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -3.231952	    1.518588
[ Info:   1	     -2.838866	    0.275166
[ Info:   2	     -2.843659	    0.048131
[ Info:   3	     -2.835333	    0.009764
[ Info:   4	     -2.835105	    0.001849
[ Info:   5	     -2.834863	    0.000380
[ Info:   6	     -2.834844	    0.000076
[ Info:   7	     -2.834837	    0.000016
[ Info:   8	     -2.834836	    0.000003
[ Info:   9	     -2.834836	    0.000001
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
using Slater Transition State approach.
```julia-repl
julia> const ev = 27.21138400
27.211384

julia> EA = -lda(Z = 3, conf = ((2, 1.5),()), β = 0.8).orbitals[3,end]*ev
[ Info: Using logarithmic 401 point grid and step dx = 0.1250
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000 
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -6.539026	    3.615572
[ Info:   1	     -7.959516	    2.552225
[ Info:   2	     -7.352136	    0.795066
[ Info:   3	     -7.361211	    0.204677
[ Info:   4	     -7.367748	    0.036300
[ Info:   5	     -7.367463	    0.007819
[ Info:   6	     -7.367563	    0.001583
[ Info:   7	     -7.367568	    0.000356
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
0.7269491136274684

julia> IP = -lda(Z = 3, conf = ((2, 0.5),()), β = 0.8).orbitals[3,end]*ev
[ Info: Using logarithmic 401 point grid and step dx = 0.1250
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000 
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle		energy		|Δρ|
[ Info:   0	     -8.080241	    1.624665
[ Info:   1	     -7.297013	    0.316713
[ Info:   2	     -7.261408	    0.082944
[ Info:   3	     -7.259677	    0.021616
[ Info:   4	     -7.259888	    0.005624
[ Info:   5	     -7.260061	    0.001467
[ Info:   6	     -7.260127	    0.000381
[ Info:   7	     -7.260148	    0.000099
[ Info:   8	     -7.260154	    0.000025
[ Info:   9	     -7.260155	    0.000007
[ Info:  10	     -7.260156	    0.000002
[ Info:  11	     -7.260156	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               7.168653
│       ELECTRON-ELECTRON              3.619145
│       EXCHANGE-CORRELATION          -1.588566 
│       ELECTRON-NUCLEAR             -16.459388
│       TOTAL ENERGY                  -7.260156
└       VIRIAL RATIO                   2.012764
5.316382271361109
```
Experimental values are: IP = 5.39 eV, EA = 0.62 eV.

## Hooke atom
`https://en.wikipedia.org/wiki/Hooke's_atom`

An exact solution for Hooke's atom is E = 2.0 a.u.

```julia-repl
julia> r, n, dx = radial_grid(301);

julia> lda((r, n, dx), conf = 2, vp = 1/8 * r.^2, β = 0.8).energy;
[ Info: Neutral atom with Z =  2 is assumed.
[ Info: Using logarithmic 301 point grid and step dx = 0.1667
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000 
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle		energy		|Δρ|
[ Info:   0	      2.103849	    1.896284
[ Info:   1	      2.013080	    0.358056
[ Info:   2	      2.022087	    0.074260
[ Info:   3	      2.025397	    0.014855
[ Info:   4	      2.026076	    0.003008
[ Info:   5	      2.026201	    0.000622
[ Info:   6	      2.026224	    0.000131
[ Info:   7	      2.026228	    0.000028
[ Info:   8	      2.026229	    0.000006
[ Info:   9	      2.026229	    0.000001
[ Info:  10	      2.026229	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022578
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
function lda(grid = radial_grid();
                Z = nothing, 
             conf = nothing,
               xc = SVWN,
               vp = nothing,
             ρ_in = nothing,
                β = 0.3,
               δn = 1.0e-6,
               δE = 5.0e-7,
            maxit = 50)
  
  if Z == nothing && conf == nothing 
    @error "System is not specified: provide ether atomic number Z 
       or/and electron configuration.
       
       Examples:
       --------
       
       # Helium
       lda(Z = 2)
       # Neon
       lda(conf = atomic_electron_configuration[:Ne])
       # Li+
       lda(Z = 3, conf = atomic_electron_configuration[:He])
       # He 1s1 2s1 excited state
       lda(Z = 2, conf = ((1, 1), ())"
    return (energy = nothing, density = nothing, orbitals = nothing)
  
  elseif Z == nothing
    Q = sum(Iterators.flatten(conf))
    Z  = Q  
    @info @sprintf("Neutral atom with Z = %2i is assumed.", Z)
  
  elseif conf == nothing
    conf = atomic_electron_configuration[Z]
    Q = Z
    s = keys(atomic_electron_configuration)[Z]
    @info @sprintf("Neutral atom with Z = %2i and %s electron configuration is assumed.", Z, s)
  else
    Q = sum(Iterators.flatten(conf))
  end
  
  r, n, dx = grid
  @info @sprintf("Using logarithmic %3i point grid and step dx = %5.4f", n, dx)
  
  if vp == nothing vp = -Z ./ r end

  if ρ_in == nothing   
    ρ_in = TF.(r, Q)
    @info "Using Thomas-Fermi starting electron density"
  end
  
  r² = r .* r; r³ = r² .* r; sqr = sqrt.(r)
  
  nlevels = length(collect(Iterators.flatten(conf)))
  ψ = zeros(n + 4, nlevels)

  Δ = laplacian(n, dx)
  H = -1/2*Δ; for i=1:n H[i,:] ./= r²[i] end
  T = diag(H)
  Δ[diagind(Δ, 0)] .-= 1/4;

  ev_sum = 0.0; E = 0.0; Δρ = 0.0
  ρ_out, vh, vxc, εxc = [Array{Float64,1}(undef, n) for _ = 1:4]

  @info @sprintf("Starting SCF procedure with density mixing parameter β = %5.4f 
      and convergence threshold |Δρ| ≤ %e", β, δn)
  @info "cycle\t\tenergy\t\t|Δρ|"
  for i = 0:maxit
    # solve Laplace equation to find Hartree potential
    
    vh = Δ \ (-4π * ρ_in .* sqr .* r²)
   
    # apply boundary conditions for vh at r → 0 and r → ∞
   
    vh .-= (vh[end] - Q / sqr[end])/sqr[end] .* sqr .+
           (vh[  1] - Q * sqr[  1])*sqr[  1] ./ sqr
    vh ./= sqr
   
    # calculate DFT exchange-correlation potential vxc and energy density εxc
    # see https://www.theoretical-physics.net/dev/quantum/dft.html#the-xc-term
   
    for i=1:n vxc[i], εxc[i] = xc(ρ_in[i]) end
   
    # assemble Kohn-Sham potential
    # see https://www.theoretical-physics.net/dev/quantum/dft.html#kohn-sham-equations
   
    V = vp + vh + vxc
   
    # solve Schrödinger equation to find new density - ρ_out
    # and one-particle energy ev_sum
   
    H[diagind(H, 0)] = T
    ev_sum = radial_shr_eq!(grid, V, conf, H, ρ_out, ψ)
   
    # DFT total energy:
    # see https://www.theoretical-physics.net/dev/quantum/dft.html#total-energy
    
    E_prev = E
    E = ev_sum + 4π * ∫(dx, ρ_out .* (-1/2 * vh .- vxc .+ εxc) .* r³)
      
    # no more charge oscillations - density converged
   
    Δρ = 4π * ∫(dx, abs.(ρ_out - ρ_in) .* r³)

    if Δρ < δn && abs(E_prev - E) < δE
      @info @sprintf "%3i\t%14.6f\t%12.6f\n" i E Δρ
      break
    end  
   
    # admix density from previous iteration to converge SCF
    
    ρ_in = (1.0 - β) * ρ_in + β * ρ_out
    
    @info @sprintf "%3i\t%14.6f\t%12.6f\n" i E Δρ
  end
  Δρ > δn && @warn "SCF did not converge after $maxit iterations!"
 
  # analysis of energy components
  denoise!(ρ_out)

  Ek  = ev_sum - 4π * ∫(dx, ρ_out .* ( vp + vh + vxc ) .* r³)
  Eh  = 2π * ∫(dx, ρ_out .* vh .* r³)
  Exc = 4π * ∫(dx, ρ_out .* εxc .* r³)
  Ep  = 4π * ∫(dx, ρ_out .* vp .* r³)
  VR  = (Eh + Exc + Ep) / -Ek

  @info @sprintf("RESULTS SUMMARY:
      ELECTRON KINETIC       %16.6f
      ELECTRON-ELECTRON      %16.6f
      EXCHANGE-CORRELATION   %16.6f 
      ELECTRON-NUCLEAR       %16.6f
      TOTAL ENERGY           %16.6f
      VIRIAL RATIO           %16.6f\n", Ek, Eh, Exc, Ep, E, VR) 
  return (energy = E, density = ρ_out, orbitals = ψ)
end

end # module
