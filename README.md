# AtomEnergyLevels

Solve numerically 

 1. 1D [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) by the [spectral collocation (i.e., pseudospectral) method](https://en.wikipedia.org/wiki/Collocation_method)
 2. Schrödinger equation for the [particle in a spherically symmetric potential](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)
 3. [Kohn-Sham equation](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations) for an atom using [local-density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

## Usage

### Morse potential

Using the following code, one can find the vibrational levels of the radical OH⋅, 
for which the potential energy surface is approximated through the Morse potential.

The bulk of the calculations were performed by the Laplacian function.

```
  julia> using LinearAlgebra, AtomEnergyLevels  
  julia> n = 501;  
  julia> D = 0.1994; β = 1.189; x₀ = 1.821;  
  julia> mH = 1.00794; mO = 15.9994; μ = 1822.8885*(mH*mO)/(mH+mO);  
  julia> xmin = 2.0 - x₀; xmax = 12.0 + x₀;  
  julia> x = range(xmin, xmax, length = n); dx = x.step.hi;  
  julia> V(x) = D*(exp(-β*(x-x₀))-1)^2 - D;  # Morse potential  
  julia> Δ = laplacian(n, dx);  
  julia> H = -1/2μ*Δ + Diagonal(V.(x));      # Hamiltonian  
  julia> ϵ, ψ = eigen(H);                    # numerical solution  
  julia> ω = β*sqrt(2D/μ); δ = ω^2 / 4D;  
  julia> E(n) = ω*(n+1/2) - δ*(n+1/2)^2 - D; # exact solution
  julia> ϵ[1]
  -0.19047201661032384
  julia> E(0)
  -0.19047201661032412  
```

### 3D isotropic harmonic oscillator problem

Using the following code, one can find the energy levels for a
[spherically-symmetric three-dimensional harmonic oscillator.](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator)

To find the energy levels with quantum numbers nᵣ = 1, 2, 3 and l = 0, 1, 2, 
the levels of interest should be included in the tuple of the 
electronic configuration as follows

```
#          l=0        l=1        l=2
conf = ((0, 0, 0), (0, 0, 0), (0, 0, 0));
#   nᵣ = 0, 1, 2;   0, 1, 2;   0, 1, 2;  
```

where the tuples `(0,0,0)` indicate population numbers of each shell.

```
julia> r, n, dx = radial_grid(300, -20, 5);
julia> ψ = radial_shr_eq((r, n, dx), 1/2*r.^2, conf).orbitals;
```

The output should be a matrix

```
julia> ψ[1:3,:]
3×9 Array{Float64,2}:
 0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0  2.0  # azimuthal quantum numbers, l
 0.0  1.0  2.0  0.0  1.0  2.0  0.0  1.0  2.0  # radial quantum numbers, nᵣ
 1.5  3.5  5.5  2.5  4.5  6.5  3.5  5.5  7.5  # energy levels
```

### Hooke's atom

Find the total energy and density of the [Hooke's atom](https://en.wikipedia.org/wiki/Hooke's_atom) using LDA DFT.

```
julia> using AtomEnergyLevels
julia> r, n, dx = radial_grid(300);
julia> E, ρ, ψ = lda((r, n, dx), conf = 2, vₑₓ = 1/8 * r.^2, β = 0.8);
[ Info: Neutral atom with Z =  2 is assumed.
[ Info: Using logarithmic 300 point grid and step dx = 0.1672
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000 
└       and convergence threshold |Δρ| ≤ 1.000000e-06
[ Info: cycle		energy		|Δρ|
[ Info:   0	      2.103849	    1.904091
[ Info:   1	      2.013080	    0.355807
[ Info:   2	      2.022087	    0.074240
[ Info:   3	      2.025397	    0.014783
[ Info:   4	      2.026076	    0.002977
[ Info:   5	      2.026201	    0.000618
[ Info:   6	      2.026224	    0.000130
[ Info:   7	      2.026228	    0.000027
[ Info:   8	      2.026229	    0.000006
[ Info:   9	      2.026229	    0.000001
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022578
│       EXCHANGE-CORRELATION          -0.523773 
│       ELECTRON-NUCLEAR               0.899965
│       TOTAL ENERGY                   2.026229
└       VIRIAL RATIO                  -2.229260
julia> using SpecialFunctions, PyPlot
julia> N² = 1 / (π^(3/2) * (8 + 5 * sqrt(π)));
julia> ρₑ = @. 2N² * exp(-1/2 * r^2) * (sqrt(π/2) * (7/4 + 1/4 * r^2
                         + (r + 1/r) * erf(r/sqrt(2))) + exp(-1/2 * r^2));
julia> pygui(true);
julia> plot(r, ρₑ, label = "exact"); plot(r, ρ, label = "LDA");
```

![Comparison of the exact Hooke's atom density with LDA numerical result](./hooke_atom_density.png)

## Author

Sergei Malykhin, s.e.malykhin@gmail.com

## License

This project is licensed under the MIT License - see the LICENSE file for
details.
