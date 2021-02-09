# AtomEnergyLevels

[![Build Status](https://travis-ci.com/malykhin-sergei/AtomEnergyLevels.jl.svg?branch=master)](https://travis-ci.com/malykhin-sergei/AtomEnergyLevels.jl)
[![Coverage Status](https://coveralls.io/repos/github/malykhin-sergei/AtomEnergyLevels.jl/badge.svg?branch=master)](https://coveralls.io/github/malykhin-sergei/AtomEnergyLevels.jl?branch=master)
[![codecov](https://codecov.io/gh/malykhin-sergei/AtomEnergyLevels.jl/branch/master/graph/badge.svg?token=SZG44ANDAN)](https://codecov.io/gh/malykhin-sergei/AtomEnergyLevels.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://malykhin-sergei.github.io/AtomEnergyLevels.jl/dev/)

Solve numerically

 1. 1D [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) by the [spectral collocation (i.e., pseudospectral) method](https://en.wikipedia.org/wiki/Collocation_method)
 2. Schrödinger equation for the [particle in a spherically symmetric potential](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)
 3. [Kohn-Sham equation](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations) for an atom using [local-density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

## Usage

### Morse potential

Using the following code, one can find the vibrational levels of the radical OH⋅,
for which the potential energy surface is approximated through the [Morse potential](https://en.wikipedia.org/wiki/Morse_potential).

```julia
using AtomEnergyLevels
using LinearAlgebra, Test, Printf

# Morse potential
V(x, D, β, x₀) = D*(exp(-β*(x-x₀))-1)^2 - D;

# parameters of the O-H bond
D  = 0.1994;   β = 1.189;  x₀ = 1.821;
mH = 1.00794; mO = 15.9994; μ = 1822.8885*(mH*mO)/(mH+mO);  

# grid
x = (2.0 - x₀):0.1:(12.0 + x₀);

# hamiltonian
H = -1/2μ*laplacian(x) + Diagonal(V.(x, D, β, x₀));

# numerical solution
ϵ, ψ = eigen(H);                      

# exact solution
ω = β*sqrt(2D/μ); δ = ω^2 / 4D;  
E(n) = ω*(n+1/2) - δ*(n+1/2)^2 - D;

# comparison
@testset "energies of the Morse potential" begin
  for i in 1:5
    @printf "Level %i: E(exact) = %5.10f E(approx) = %5.10f\n" i E(i-1) ϵ[i]
    @test ϵ[i] ≈ E(i-1)
  end
end
```
Output
```
Level 1: E(exact) = -0.1904720166 E(approx) = -0.1904720166
Level 2: E(exact) = -0.1732294768 E(approx) = -0.1732294768
Level 3: E(exact) = -0.1568048397 E(approx) = -0.1568048397
Level 4: E(exact) = -0.1411981052 E(approx) = -0.1411981048
Level 5: E(exact) = -0.1264092734 E(approx) = -0.1264092731
Test Summary:                   | Pass  Total
energies of the Morse potential |    5      5
Test.DefaultTestSet("energies of the Morse potential", Any[], 5, false)
```
### 3D isotropic harmonic oscillator problem

Using the following code, one can find the energy levels for a
[spherically-symmetric three-dimensional harmonic oscillator.](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator)

To find the energy levels with quantum numbers nᵣ = 0, 1, 2 and l = 0, 1, 2,
the levels of interest should be included in the electronic configuration
```julia
using AtomEnergyLevels, Printf

function isotropic_harmonic_oscillator(config)
  ψ = radial_shr_eq(r -> 1/2*r^2, conf = config).orbitals;
  @printf("\tϵ(calc.)\tϵ(exact)\tΔϵ\n");
  for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
    nᵣ, l = quantum_numbers
    ϵ_calc = orbital.ϵᵢ
    n = nᵣ + l + 1
    ϵ_exact = 2nᵣ + l + 3/2
    @printf("%i%s\t%10.8f\t%10.8f\t%+0.6e\n",
            n, shells[l+1], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)
  end
end

isotropic_harmonic_oscillator(c"1s1 2s1 3s1 2p1 3p1 4p1 3d1 4d1 5d1");
```
The output will be
```
        ϵ(calc.)        ϵ(exact)        Δϵ
1S      1.50000000      1.50000000      -2.639289e-11
2P      2.50000000      2.50000000      +3.940093e-11
3D      3.50000000      3.50000000      +2.014788e-11
2S      3.50000000      3.50000000      -2.537082e-12
3P      4.50000000      4.50000000      +3.407141e-11
4D      5.50000000      5.50000000      +3.625988e-11
3S      5.50000000      5.50000000      +4.071410e-12
4P      6.50000000      6.50000000      +2.108536e-12
5D      7.50000000      7.50000000      +4.278622e-11
```
### Hooke's atom

Find the total energy and density of the [Hooke's atom](https://en.wikipedia.org/wiki/Hooke's_atom) using LDA DFT.

```julia
using AtomEnergyLevels
using SpecialFunctions, PyPlot
pygui(true);

x = -30.0:0.1:20.0;
r  = exp.(x);

# exact density of the ground state, see
# S. Kais et al // Density functionals and dimensional 
# renormalization for an exactly solvable model, JCP 99, 417 (1993); 
# http://dx.doi.org/10.1063/1.465765

N² = 1 / (π^(3/2) * (8 + 5 * sqrt(π)));
ρₑ = @. 2N² * exp(-1/2 * r^2) * (sqrt(π/2) * (7/4 + 1/4 * r^2 + (r + 1/r) * erf(r/sqrt(2))) + exp(-1/2 * r^2));

Xα(alpha = 1) = xc_func(xc_type_lda, false, ρ -> LDA_X(ρ, α = alpha), ρ -> 0)

xalpha = lda(2, x, conf = c"[He]", Vex = r -> 1/8 * r^2, xc = Xα(0.798));
svwn5  = lda(2, x, conf = c"[He]", Vex = r -> 1/8 * r^2);

begin
    title("Hooke's atom density")
    xlabel(L"\frac{r Z^{1/3}}{b},\,a.u.")
    ylabel(L"\rho Z^{-\frac{4}{3}}")
    
    ax1 = PyPlot.axes()
    ax1.set_xlim([0,5])

    plot(r, ρₑ, label = "Exact")
    plot(r, xalpha.density ./ r, label = "Xα, α=0.798")
    plot(r, svwn5.density ./ r, label = "LDA")

    legend(loc = "upper right", fancybox = "true")
    grid("on")
end
```
LDA results output:
```
[ Info: A grid of 501 points x = -30.0000:0.1000:20.0000 is used
[ Info: Starting SCF procedure with convergence threshold δn ≤ 1.490116e-08
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
```
Compare with exact Kohn-Sham values from Table II of S. Kais et al
```
        ELECTRON KINETIC               0.6352
        ELECTRON-ELECTRON              1.0320
        EXCHANGE-CORRELATION          -0.5553
        ELECTRON-NUCLEAR               0.8881
        TOTAL ENERGY                   2.0000
        VIRIAL RATIO                  -2.1486
```
![Comparison of the exact Hooke's atom density with LDA numerical result](./examples/hooke_atom.png)

### Uranium atom

Compare with NIST [Atomic Reference Data for Electronic Structure
Calculations, Uranium](https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-90)

```julia
using AtomEnergyLevels
using PyPlot
pygui(true)

b = 1/2 * (3π/4)^(2/3)

x = -30.0:0.1:20.0; r  = exp.(x);

Z = 92; ρ = lda(Z, x).density;

begin
    title("Radial electronic density of the U and Thomas-Fermi atoms")
    xlabel(L"\frac{r Z^{1/3}}{b},\,a.u.")
    ylabel(L"\rho Z^{-\frac{4}{3}}")

    grid("on")

    ax1 = PyPlot.axes()
    ax1.set_xlim([0.0001,100])
    PyPlot.xscale("log")
    
    plot(r/b, 4π * r .^ 2 .* TF.(r, 1), label = "TF")
    plot(r/b * Z^(1/3), 4π * r .* ρ / Z^(4/3), label = "U (LDA)")
    legend(loc = "upper right", fancybox = "true")
end
```
Output:
```
[ Info: A grid of 501 points x = -30.0000:0.1000:20.0000 is used
[ Info: Starting SCF procedure with convergence threshold δn ≤ 1.490116e-08
[ Info: cycle           energy          δn
[ Info:   0      -25892.940688     13.94682700
[ Info:   1      -25392.907488      9.14239485
[ Info:   2      -25958.273128     10.84751085
[ Info:   3      -25323.130397     13.21089613
[ Info:   4      -26013.215074     13.11048291
[ Info:   5      -25360.517609     11.72644312
[ Info:   6      -25650.759429      4.29309894
[ Info:   7      -25662.547556      1.93857182
[ Info:   8      -25660.046452      0.76699314
[ Info:   9      -25658.557814      0.24604636
[ Info:  10      -25658.398334      0.03850339
[ Info:  11      -25658.437820      0.00531584
[ Info:  12      -25658.394604      0.00138062
[ Info:  13      -25658.399835      0.00115990
[ Info:  14      -25658.417355      0.00041016
[ Info:  15      -25658.417852      0.00023405
[ Info:  16      -25658.417853      0.00021146
[ Info:  17      -25658.417854      0.00019107
[ Info:  18      -25658.417855      0.00017267
[ Info:  19      -25658.417857      0.00015605
[ Info:  20      -25658.417858      0.00014105
[ Info:  21      -25658.417859      0.00012749
[ Info:  22      -25658.417861      0.00011525
[ Info:  23      -25658.417874      0.00001674
[ Info:  24      -25658.417887      0.00000285
[ Info:  25      -25658.417890      0.00000071
[ Info:  26      -25658.417887      0.00000031
[ Info:  27      -25658.417889      0.00000005
[ Info:  28      -25658.417888      0.00000003
[ Info:  29      -25658.417889      0.00000001
[ Info: SCF has converged after 29 iterations
┌ Info: RESULTS SUMMARY:
│       ELECTRON KINETIC           25651.231179
│       ELECTRON-ELECTRON           9991.594177
│       EXCHANGE-CORRELATION        -425.032628
│       ELECTRON-NUCLEAR          -60876.210617
│       TOTAL ENERGY              -25658.417889
└       VIRIAL RATIO                   2.000280
```
![Radial electronic density of the U and Thomas-Fermi atoms](./examples/shell_structure.png)

## Author

Sergei Malykhin, s.e.malykhin@gmail.com

## License

This project is licensed under the MIT License - see the LICENSE file for
details.