# AtomEnergyLevels.jl

Solve numerically

 1. 1D [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) by the [spectral collocation (i.e., pseudospectral) method](https://en.wikipedia.org/wiki/Collocation_method)
 2. Schrödinger equation for the [particle in a spherically symmetric potential](https://en.wikipedia.org/wiki/Particle_in_a_spherically_symmetric_potential)
 3. [Kohn-Sham equation](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations) for an atom using [local-density approximation (LDA)](https://en.wikipedia.org/wiki/Local-density_approximation)

## Theory

### Change of variables

The dynamics of a particle in a spherically symmetric potential are governed
by a Hamiltonian of the following form:
```math
-\frac{\hbar^{2}}{2\mu}\underset{\Delta R_{nl}}{\underbrace{\left[\frac{1}{r^{2}}\frac{\partial}{\partial r}\left(r^{2}\frac{\partial R_{nl}}{\partial r}\right)-\frac{l(l+1)}{r^{2}}R_{nl}(r)\right]}}+V(r)R_{nl}(r)=E_{nl}R_{nl}(r)
```
where ``\Delta`` - [Laplace operator](https://en.wikipedia.org/wiki/Laplace_operator),
``r`` - distance between the particle and a defined center point (nucleus),
``\mu`` - [reduced mass](https://en.wikipedia.org/wiki/Reduced_mass),
``R_{nl}(r)`` - radial part of a wavefunction,
``E_{nl}`` - eigenstate of a system (energy level).
The eigenfunctions take the form
```math
\psi(r,\,\theta,\,\phi)=R_{nl}(r)\cdot Y_{lm}(\theta,\,\phi)
```
where ``Y_{lm}(\theta,\,\phi)`` - [spherical harmonic.](https://en.wikipedia.org/wiki/Spherical_harmonics)

The substitution
```math
R_{nl}(r)=\frac{1}{r}P_{nl}(r)
```
results in Schrödinger equation of the form:
```math
-\frac{\hbar^{2}}{2m}\left[\frac{\partial^{2}P_{nl}}{\partial r^{2}}-\frac{l(l+1)}{r^{2}}P_{nl}(r)\right]+V(r)P_{nl}(r)=E_{nl}P_{nl}(r)
```
To solve the latter equation the variable ``r`` and eigenfunction ``P_{nl}(r)``
are changed to ``x`` and ``y_{nl}(x)``, where

```math
\begin{aligned}
r &= e^{x}\\
y_{nl}(x) &= \frac{1}{\sqrt{r}}P_{nl}\left(r(x)\right)
\end{aligned}
```

After that, radial Schrödinger equation takes convenient form for numerical solution.
```math
-\frac{1}{2\mu}\frac{\partial^{2}y_{nl}(x)}{\partial x^{2}}+\left(\frac{1}{2\mu}\left(l+\frac{1}{2}\right)^{2}+r^{2}V(r)\right)y_{nl}(x)=r^{2}E_{nl}y_{nl}(x)
```

[Atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used
throughout, i.e. ``\hbar = 1`` and ``4\pi\epsilon_0 = 1``. It means that
distances are in 1 a.u. = 0.529177 Angstrom and coulomb potential is
```math
V(r) = -\frac{Z}{\epsilon r}
```

### Pseudospectral method

[Pseudospectral method](https://en.wikipedia.org/wiki/Collocation_method) is
convenient approach to code a solution of differential equation. The wavefunction
is taken as a linear combination of [sinc-functions.](https://en.wikipedia.org/wiki/Sinc_function)
```math
\begin{aligned}
y(x) &= \sum_{i=1}^{n}c_{i}\phi_{i}(x)\\
\phi_{i}(x) &= \frac{\sin\left(\pi(x-x_{i})/h\right)}{\pi(x-x_{i})/h}
\end{aligned}
```
Discrete grid
```math
x_{i}=x_{min}+i\cdot h
```
where ``i=1\ldots n`` and step size ``h``, defines basis set, since
```math
\phi_{i}(x_{j})=\left\{ \begin{array}{clr}
1 & \text{for} & i=j\\
0 & \text{for} & i\neq j
\end{array}\right.
```
Coefficients are values of the function on a grid.
```math
c_i = y(x_i)
```
The first derivative of a basis function taken on a grid point is
```math
\phi'_j(x_i)=\frac{1}{h}\frac{(-1)^i}{i}
```
Therefore, differentiation of a function on the grid ``x_i`` is
equivalent to matrix-vector product
```math
\mathrm{y'(x_i)} = \mathrm{D^{(1)}} \cdot \mathrm{y(x_i)}
```
where matrices for the first and second derivatives are
```math
\mathrm{D^{(1)}}=\frac{1}{h}\left(\begin{array}{ccccc}
0 & 1 & -\frac{1}{2} & \cdots & \frac{(-1)^{n}}{n-1}\\
-1 & 0 & 1\\
\frac{1}{2} & -1 & 0 &  & -\frac{1}{2}\\
 &  &  & \ddots & 1\\
\frac{(-1)^{n-1}}{n-1} & \cdots & \frac{1}{2} & -1 & 0
\end{array}\right)
```
```math
\mathrm{D^{(2)}}=\frac{1}{h^{2}}\left(\begin{array}{ccccc}
-\frac{\pi^{2}}{3} & 2 & -\frac{1}{2} & \cdots & \frac{2(-1)^{n}}{(n-1)^{2}}\\
2 & -\frac{\pi^{2}}{3} & 2\\
-\frac{1}{2} & 2 & -\frac{\pi^{2}}{3} &  & -\frac{1}{2}\\
 &  &  & \ddots & 2\\
\frac{2(-1)^{n}}{(n-1)^{2}} & \cdots & -\frac{1}{2} & 2 & -\frac{\pi^{2}}{3}
\end{array}\right)
```
With these matrices, the numerical solution of Schrödinger equation becomes
a linear algebra problem. For this we have the following ingredients:

1. logarithmic grid ``r_i=\exp(x_i)`` on which vectors are calculated: orbitals ``y_k`` and potential ``V(r)``,
2. diagonal matrices: ``\mathrm{S} = r_i^2`` and ``\mathrm{V_{\text{eff}}}=\frac{1}{2\mu}\left(l+\frac{1}{2}\right)^{2}+r_i^{2}V(r_i)``,  
3. kinetic energy operator ``\mathrm{K}=-\frac{1}{2\mu}\cdot D^{(2)}``,
4. Hamiltonian operator ``\mathrm{H = K + V_{\text{eff}}}``.

Variants of this approach can be found in the documentation
[A Matlab Differentiation Matrix Suite.](http://appliedmaths.sun.ac.za/~weideman/research/differ.html)

### Peculiarities of the eigenvalue problem

Solution of the [eigenvalues problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)
is pairs of eigenvectors ``\epsilon_k`` and eigenvectors ``y_k``,
that correspond to the particle's energy levels and wavefunctions (orbitals).
The ``i`` component of the ``y_k`` eigenvector is the value of the ``k`` level
wavefunction at the ``x_i`` grid point. It is this property of the
pseudospectral method that makes the algorithm compact and clear.
```math
\mathrm{H\cdot y = S\cdot y}
```
The [generalized eigenvalue problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem)
of eigenvalues of a symmetrical matrix ``H`` and a symmetrical,
positively defined matrix ``S`` is reduced to the symmetrical eigenvalues
problem using the [symmetrical Löwdin orthogonalization.](https://doi.org/10.1016/S0065-3276(08)60339-1)
```math
\mathrm{\underset{H'}{\underbrace{S^{-\frac{1}{2}}\cdot H\cdot S^{-\frac{1}{2}}}}\cdot\underset{y'}{\underbrace{S^{\frac{1}{2}}\cdot y}}=E\cdot\underset{y'}{\underbrace{S^{\frac{1}{2}}\cdot y}}}
```
where
```math
\mathrm{H'}=\left(\begin{array}{cccc}
\frac{H_{11}}{r_{1}r_{1}} & \frac{H_{12}}{r_{1}r_{2}} & \cdots & \frac{H_{1n}}{r_{1}r_{n}}\\
\frac{H_{21}}{r_{2}r_{1}} & \frac{H_{22}}{r_{2}^{2}} &  & \vdots\\
\vdots &  & \ddots\\
\frac{H_{n1}}{r_{n}r_{1}} & \cdots &  & \frac{H_{nn}}{r_{n}^{2}}
\end{array}\right)
```
However, the floating point representation of numbers leads to a loss of
symmetry of the matrix ``H'`` due to non-associativity.
```math
H_{ij}/r_i/r_j \neq H_{ji}/r_j/r_i
```
[The standard algorithm implemented in the LAPACK library](http://www.netlib.org/lapack/lug/node54.html)
for the generalized eigenvalue problem uses the [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition),
but its scope [is limited](http://www.netlib.org/utk/people/JackDongarra/etemplates/node177.html#sec:gsym_pert)
to the case when the ``\mathrm{S}`` matrix is well defined, i.e.
condition number
```math
k(\mathrm{S})\equiv\left\Vert \mathrm{S}\right\Vert _{2}\cdot\left\Vert \mathrm{S}^{-1}\right\Vert _{2}
```
is not large. In our case, this number is very large, ca ``10^{52}``
```math
k(\mathrm{S})=\sqrt{\sum_{i}r_{i}^{2}}\cdot\sqrt{\sum_{i}\frac{1}{r_{i}^{2}}}>\frac{r_{max}}{r_{min}}
```
To remedy the problem, note that eigenvalues of [matrix pencils](https://en.wikipedia.org/wiki/Matrix_pencil)
``(\mathrm{H,S})`` and ``(\mathrm{H,\alpha S + \beta H})`` are same, when
matrix ``\mathrm{\alpha S + \beta H}`` is positive definite. Then
```math
\mathrm{H \cdot y = \theta \cdot (\alpha S + \beta H) \cdot y}
```
has eigenvalues
```math
\theta_i = \frac{\epsilon_i}{\alpha + \beta \epsilon_i}
```
Setting ``\beta = 1`` and ``\alpha = 10^5`` we make new matrix
```math
\mathrm{S' = \alpha S + H}
```
which is positive definite and well-conditioned, then solve
```math
\mathrm{H \cdot y = \theta \cdot S' \cdot y}
```
using standard LAPACK routine. Eigenvalues of interest are found with
```math
\epsilon_i = \frac{\alpha \cdot \theta_i}{1 - \theta_i}
```
This approach is used here.

### Sketch of the solution

We are going to find a solution of the differential equation
```math
-\frac{1}{2\mu}\frac{\partial^{2}y_{nl}(x)}{\partial x^{2}}+\left(\frac{1}{2\mu}\left(l+\frac{1}{2}\right)^{2}+r^{2}V(r)\right)y_{nl}(x)=r^{2}E_{nl}y_{nl}(x)
```
on the grid ``x_i = x_{min} + (i-1)\cdot dx``.

```@example
using LinearAlgebra

# find differentiation matrix
function laplacian(x)
  n, h = length(x), step(x)
  Δ = Matrix{Float64}(undef, n, n)
  Δ[diagind(Δ, 0)] .= -1/3 * π^2 / h^2
  for i = 2:n
    Δ[diagind(Δ, i - 1)] =
    Δ[diagind(Δ, 1 - i)] .= 2 * (-1)^i / (i - 1)^2 / h^2
  end
  return Δ
end

# make logarithmic grid on r
x = -30.0:0.175:5.0; r = exp.(x); r² = exp.(2x)

# set up parameters of the hydrogen atom hamiltonian:
# reduced mass, azimuthal quantum number, and nucleus charge
μ = 1; l = 0; Z = 1;

# make hamiltonian matrix
H = -1/2μ*laplacian(x) + Diagonal(1/2μ * (l + 1/2)^2 .- Z*r);

# form positive definite well-conditioned matrix S
α = 1e5; β = 1;
S = β*H + α * Diagonal(r .* r);

# solve generalized eigenproblem
θ, y = eigen(H, S);
ϵ = α * θ ./ (1.0 .- θ);

# see energy of the 1s-state of hydrogen atom
ϵ[1]
```
## Functions for solving Schrödinger equation
### Differentiation matrix

```@docs
laplacian
```

Using the following code, one can find the vibrational levels of the radical OH⋅,
for which the potential energy surface is approximated through the [Morse potential](https://en.wikipedia.org/wiki/Morse_potential).

```@example
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
for i in 1:5
  @printf "Level %i: E(exact) = %5.10f E(approx) = %5.10f\n" i E(i-1) ϵ[i]
end
```
### Atomic electron configuration parsing
```@docs
conf_enc
```
### Radial Schrödinger equation solver
```@docs
radial_shr_eq
```

#### Isotropic harmonic oscillator

Using the following code, one can find the energy levels for a
[spherically-symmetric three-dimensional harmonic oscillator.](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Example:_3D_isotropic_harmonic_oscillator)

To find the energy levels with quantum numbers nᵣ = 0, 1, 2 and l = 0, 1, 2
the levels of interest should be included in the electronic configuration.

```@example
using AtomEnergyLevels, Printf

function isotropic_harmonic_oscillator(cfg)
  ψ = radial_shr_eq(r -> 1/2*r^2, conf = conf_enc(cfg)).orbitals;
  @printf("\tϵ(calc.)\tϵ(exact)\tΔϵ\n");
  for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
    nᵣ, l = quantum_numbers
    ϵ_calc = orbital.ϵᵢ
    n = nᵣ + l + 1
    ϵ_exact = 2nᵣ + l + 3/2
    @printf("%i%s\t%10.8f\t%10.8f\t%+0.6e\n",
            n, atomic_shell[l], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)
  end
end

isotropic_harmonic_oscillator("1s1 2s1 3s1 2p1 3p1 4p1 3d1 4d1 5d1");
```

#### Kratzer Potential

```@example
using AtomEnergyLevels, Printf

"""
Solutions of the Schrödinger Equation for the Kratzer Potential.
A. Kratzer, "Die ultraroten Rotationsspektren der Halogenwasserstoffe," 
Zeitschrift für Physik, 3(5), 1920 pp. 289–307. doi:10.1007/BF01327754

See also 
https://demonstrations.wolfram.com/ExactSolutionsOfTheSchroedingerEquationForTheKratzerPotentia/
"""
function kratzer(D, a, levels::Union{UnitRange{Int64}, Int64} = 0)
    cfg = join(collect(levels) .+ 1, "s1 ") * "s1"
    ψ = radial_shr_eq(r -> -2D*(a/r - a^2 / 2r^2), 
        conf = conf_enc(cfg, maxn = last(levels) + 1)).orbitals
    
    μ = 1/2*sqrt(1 + 8a^2*D)

    @printf("\tϵ(calc.)\tϵ(exact)\tΔϵ\n");
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      ϵ_calc = orbital.ϵᵢ
      n = nᵣ + l + 1
      ϵ_exact = -2a^2 * D^2 / (nᵣ + μ + 1/2)^2
      @printf("%i%s\t%11.8f\t%11.8f\t%+0.6e\n",
              n, atomic_shell[l], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)
    end    
end

kratzer(2.5, 1.25, 0:10)
```

#### Pseudoharmonic Potential

```@example
using AtomEnergyLevels, Printf

"""
Solutions of the Schrödinger Equation for Pseudoharmonic Potential.

I.I. Gol'dman and V.D. Krivchenkov, Problems in Quantum Mechanics 
(B. T. Geǐlikman, ed., E. Marquit and E. Lepa, trans.), Reading, MA: Addison-Wesley, 1961.‬

See also
https://demonstrations.wolfram.com/ExactSolutionsOfTheSchroedingerEquationForPseudoharmonicPote/
"""
function pseudoharmonic(D, a, levels::Union{UnitRange{Int64}, Int64} = 0)
    cfg = join(collect(levels) .+ 1, "s1 ") * "s1"
    ψ = radial_shr_eq(r -> D*(r/a - a/r)^2, -5.0:0.01:4.0,
        conf = conf_enc(cfg, maxn = last(levels) + 1)).orbitals

    @printf("\tϵ(calc.)\tϵ(exact)\tΔϵ\n");
    for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers
      ϵ_calc = orbital.ϵᵢ
      n = nᵣ + l + 1
      ϵ_exact = sqrt(D/2)/a * (2 + 4nᵣ - 2a * sqrt(2D) + sqrt(1 + 8D * a^2))
      @printf("%i%s\t%11.8f\t%11.8f\t%+0.6e\n",
              n, atomic_shell[l], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)
    end    
end

pseudoharmonic(1, 2, 0:10)
```
