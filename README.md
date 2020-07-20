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
julia> r, n, dx = radial_grid(501, -30, 2);
julia> result = lda((r, n, dx), conf = 2, vp = 1/8 * r.^2, β = 0.8);
[ Info: Neutral atom with Z =  2 is assumed.
[ Info: Using logarithmic 501 point grid and step dx = 0.0640
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.8000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	      2.108127	    1.781135
[ Info:   1	      2.014134	    0.334518
[ Info:   2	      2.022496	    0.069493
[ Info:   3	      2.025678	    0.013863
[ Info:   4	      2.026331	    0.002813
[ Info:   5	      2.026450	    0.000584
[ Info:   6	      2.026473	    0.000123
[ Info:   7	      2.026477	    0.000026
[ Info:   8	      2.026478	    0.000006
[ Info:   9	      2.026478	    0.000001
[ Info:  10	      2.026478	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC               0.627459
│       ELECTRON-ELECTRON              1.022827
│       EXCHANGE-CORRELATION          -0.523773
│       ELECTRON-NUCLEAR               0.899965
│       TOTAL ENERGY                   2.026478
└       VIRIAL RATIO                  -2.229656
julia> using SpecialFunctions, PyPlot
julia> N² = 1 / (π^(3/2) * (8 + 5 * sqrt(π)));
julia> ρₑ = @. 2N² * exp(-1/2 * r^2) * (sqrt(π/2) * (7/4 + 1/4 * r^2
                         + (r + 1/r) * erf(r/sqrt(2))) + exp(-1/2 * r^2));
julia> pygui(true);
julia> plot(r, ρₑ, label = "exact"); plot(r, result.density, label = "LDA");
```

![Comparison of the exact Hooke's atom density with LDA numerical result](./hooke_atom_density.png)

### Uranium atom

```
julia> using AtomEnergyLevels, Printf
julia> @time result = lda(Z = 92, β = 0.5);
[ Info: Neutral atom with Z = 92 and U electron configuration is assumed.
[ Info: Using logarithmic 501 point grid and step dx = 0.1000
[ Info: Using Thomas-Fermi starting electron density
┌ Info: Starting SCF procedure with density mixing parameter β = 0.5000
└       and convergence threshold |Δρ| ≤ 5.000000e-07
[ Info: cycle		energy		|Δρ|
[ Info:   0	 -25892.940688	   13.946827
[ Info:   1	 -25648.785452	    6.210687
[ Info:   2	 -25663.151615	    2.715824
[ Info:   3	 -25658.522677	    1.331700
[ Info:   4	 -25658.683795	    0.663356
[ Info:   5	 -25658.531924	    0.330864
[ Info:   6	 -25658.485432	    0.165669
[ Info:   7	 -25658.455013	    0.083251
[ Info:   8	 -25658.438247	    0.041877
[ Info:   9	 -25658.428850	    0.021078
[ Info:  10	 -25658.423720	    0.010617
[ Info:  11	 -25658.420960	    0.005357
[ Info:  12	 -25658.419493	    0.002707
[ Info:  13	 -25658.418722	    0.001369
[ Info:  14	 -25658.418319	    0.000693
[ Info:  15	 -25658.418110	    0.000351
[ Info:  16	 -25658.418003	    0.000178
[ Info:  17	 -25658.417947	    0.000090
[ Info:  18	 -25658.417918	    0.000046
[ Info:  19	 -25658.417904	    0.000023
[ Info:  20	 -25658.417896	    0.000012
[ Info:  21	 -25658.417893	    0.000006
[ Info:  22	 -25658.417891	    0.000003
[ Info:  23	 -25658.417890	    0.000002
[ Info:  24	 -25658.417889	    0.000001
[ Info:  25	 -25658.417889	    0.000000
[ Info:  26	 -25658.417889	    0.000000
[ Info:  27	 -25658.417889	    0.000000
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC           25651.231179
│       ELECTRON-ELECTRON           9991.594177
│       EXCHANGE-CORRELATION        -425.032628
│       ELECTRON-NUCLEAR          -60876.210617
│       TOTAL ENERGY              -25658.417889
└       VIRIAL RATIO                   2.000280
  4.138058 seconds (29.66 k allocations: 1.552 GiB, 0.93% gc time)
julia> ψ = result.orbitals;
julia> p = sortperm(ψ[3,:]);
julia> ψ = ψ[:, p];
julia> for (i, ϵᵢ) in enumerate(ψ[3,:])
       l, nᵣ, nᵢ = ψ[1, i], ψ[2, i], ψ[4, i]
       n = nᵣ + l + 1
       @printf("\t%i%s\t(%4.1f)\t%14.6f\n", n, atomic_shell[l+1], nᵢ, ϵᵢ)
       end
	1S	( 2.0)	  -3689.355140
	2S	( 2.0)	   -639.778728
	2P	( 6.0)	   -619.108550
	3S	( 2.0)	   -161.118073
	3P	( 6.0)	   -150.978980
	3D	(10.0)	   -131.977358
	4S	( 2.0)	    -40.528084
	4P	( 6.0)	    -35.853321
	4D	(10.0)	    -27.123212
	4F	(14.0)	    -15.027460
	5S	( 2.0)	     -8.824089
	5P	( 6.0)	     -7.018092
	5D	(10.0)	     -3.866175
	6S	( 2.0)	     -1.325976
	6P	( 6.0)	     -0.822538
	5F	( 3.0)	     -0.366543
	6D	( 1.0)	     -0.143190
	7S	( 2.0)	     -0.130948
```

## Author

Sergei Malykhin, s.e.malykhin@gmail.com

## License

This project is licensed under the MIT License - see the LICENSE file for
details.
