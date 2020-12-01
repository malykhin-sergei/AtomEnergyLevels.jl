var documenterSearchIndex = {"docs":
[{"location":"#AtomEnergyLevels.jl","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Solve numerically","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"1D Schrödinger equation by the spectral collocation (i.e., pseudospectral) method\nSchrödinger equation for the particle in a spherically symmetric potential\nKohn-Sham equation for an atom using local-density approximation (LDA)","category":"page"},{"location":"#Theory","page":"AtomEnergyLevels.jl","title":"Theory","text":"","category":"section"},{"location":"#Change-of-variables","page":"AtomEnergyLevels.jl","title":"Change of variables","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"The dynamics of a particle in a spherically symmetric potential are governed  by a Hamiltonian of the following form:","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"-frachbar^22muundersetDelta R_nlunderbraceleftfrac1r^2fracpartialpartial rleft(r^2fracpartial R_nlpartial rright)-fracl(l+1)r^2R_nl(r)right+V(r)R_nl(r)=E_nlR_nl(r)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"where Delta - Laplace operator, r - distance between the particle and a defined center point (nucleus), mu - reduced mass, R_nl(r) - radial part of a wavefunction,  E_nl - eigenstate of a system (energy level). The eigenfunctions take the form","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"psi(rthetaphi)=R_nl(r)cdot Y_lm(thetaphi)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"where Y_lm(thetaphi) - spherical harmonic.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"The substitution","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"R_nl(r)=frac1rP_nl(r)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"results in Schrödinger equation of the form:","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"-frachbar^22mleftfracpartial^2P_nlpartial r^2-fracl(l+1)r^2P_nl(r)right+V(r)P_nl(r)=E_nlP_nl(r)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"To solve the latter equation the variable r and eigenfunction P_nl(r) are changed to x and y_nl(x), where","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"beginaligned\nr = e^x\ny_nl(x) = frac1sqrtrP_nlleft(r(x)right)\nendaligned","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"After that, radial Schrödinger equation takes convenient form for numerical solution.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"-frac12mufracpartial^2y_nl(x)partial x^2+left(frac12muleft(l+frac12right)^2+r^2V(r)right)y_nl(x)=r^2E_nly_nl(x)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Atomic units are used  throughout, i.e. hbar = 1 and 4piepsilon_0 = 1. It means that distances are in 1 a.u. = 0.529177 Angstrom and coulomb potential is","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"V(r) = -fracZepsilon r","category":"page"},{"location":"#Pseudospectral-method","page":"AtomEnergyLevels.jl","title":"Pseudospectral method","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Pseudospectral method is convenient approach to code a solution of differential equation. The wavefunction is taken as a linear combination of sinc-functions.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"beginaligned\ny(x) = sum_i=1^nc_iphi_i(x)\nphi_i(x) = fracsinleft(pi(x-x_i)hright)pi(x-x_i)h\nendaligned","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Discrete grid","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"x_i=x_min+icdot h","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"where i=1ldots n and step size h, defines basis set, since","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"phi_i(x_j)=left beginarrayclr\n1  textfor  i=j\n0  textfor  ineq j\nendarrayright","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Coefficients are values of the function on a grid.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"c_i = y(x_i)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"The first derivative of a basis function taken on a grid point is","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"phi_j(x_i)=frac1hfrac(-1)^ii","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Therefore, differentiation of a function on the grid x_i is equivalent to matrix-vector product","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmy(x_i) = mathrmD^(1) cdot mathrmy(x_i)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"where matrices for the first and second derivatives are","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmD^(1)=frac1hleft(beginarrayccccc\n0  1  -frac12  cdots  frac(-1)^nn-1\n-1  0  1\nfrac12  -1  0    -frac12\n      ddots  1\nfrac(-1)^n-1n-1  cdots  frac12  -1  0\nendarrayright)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmD^(2)=frac1h^2left(beginarrayccccc\n-fracpi^23  2  -frac12  cdots  frac2(-1)^n(n-1)^2\n2  -fracpi^23  2\n-frac12  2  -fracpi^23    -frac12\n      ddots  2\nfrac2(-1)^n(n-1)^2  cdots  -frac12  2  -fracpi^23\nendarrayright)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"With these matrices, the numerical solution of Schrödinger equation becomes  a linear algebra problem. For this we have the following ingredients:","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"logarithmic grid r_i=exp(x_i) on which vectors are calculated: orbitals y_k and potential V(r),\ndiagonal matrices: mathrmS = r_i^2 and mathrmV_texteff=frac12muleft(l+frac12right)^2+r_i^2V(r_i),  \nkinetic energy operator mathrmK=-frac12mucdot D^(2),\nHamiltonian operator mathrmH = K + V_texteff.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Variants of this approach can be found in the documentation  A Matlab Differentiation Matrix Suite.","category":"page"},{"location":"#Peculiarities-of-the-eigenvalue-problem","page":"AtomEnergyLevels.jl","title":"Peculiarities of the eigenvalue problem","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Solution of the eigenvalues problem  is pairs of eigenvectors epsilon_k and eigenvectors y_k,  that correspond to the particle's energy levels and wavefunctions (orbitals).  The i component of the y_k eigenvector is the value of the k level  wavefunction at the x_i grid point. It is this property of the  pseudospectral method that makes the algorithm compact and clear.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmHcdot y = Scdot y","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"The generalized eigenvalue problem of eigenvalues of a symmetrical matrix H and a symmetrical,  positively defined matrix S is reduced to the symmetrical eigenvalues  problem using the symmetrical Löwdin orthogonalization.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmundersetHunderbraceS^-frac12cdot Hcdot S^-frac12cdotundersetyunderbraceS^frac12cdot y=EcdotundersetyunderbraceS^frac12cdot y","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"where ","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmH=left(beginarraycccc\nfracH_11r_1r_1  fracH_12r_1r_2  cdots  fracH_1nr_1r_n\nfracH_21r_2r_1  fracH_22r_2^2    vdots\nvdots    ddots\nfracH_n1r_nr_1  cdots    fracH_nnr_n^2\nendarrayright)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"However, the floating point representation of numbers leads to a loss of  symmetry of the matrix H due to non-associativity.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"H_ijr_ir_j neq H_jir_jr_i","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"The standard algorithm implemented in the LAPACK library  for the generalized eigenvalue problem uses the Cholesky decomposition,  but its scope is limited  to the case when the mathrmS matrix is well defined, i.e.  condition number","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"k(mathrmS)equivleftVert mathrmSrightVert _2cdotleftVert mathrmS^-1rightVert _2","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"is not large. In our case, this number is very large, ca 10^52","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"k(mathrmS)=sqrtsum_ir_i^2cdotsqrtsum_ifrac1r_i^2fracr_maxr_min","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"To remedy the problem, note that eigenvalues of matrix pencils (mathrmHS) and (mathrmHalpha S + beta H) are same, when matrix mathrmalpha S + beta H is positive definite. Then","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmH cdot y = theta cdot (alpha S + beta H) cdot y","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"has eigenvalues","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"theta_i = fracepsilon_ialpha + beta epsilon_i","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Setting beta = 1 and alpha = 10^5 we make new matrix","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmS = alpha S + H","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"which is positive definite and well-conditioned, then solve","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"mathrmH cdot y = theta cdot S cdot y","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"using standard LAPACK routine. Eigenvalues of interest are found with","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"epsilon_i = fracalpha cdot theta_i1 - theta_i","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"This approach is used here.","category":"page"},{"location":"#Sketch-of-the-solution","page":"AtomEnergyLevels.jl","title":"Sketch of the solution","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"We are going to find a solution of the differential equation","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"-frac12mufracpartial^2y_nl(x)partial x^2+left(frac12muleft(l+frac12right)^2+r^2V(r)right)y_nl(x)=r^2E_nly_nl(x)","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"on the grid x_i = x_min + (i-1)cdot dx. ","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"using LinearAlgebra\n\n# find differentiation matrix\nfunction laplacian(x)\n  n, h = length(x), step(x)\n  Δ = Matrix{Float64}(undef, n, n)\n  Δ[diagind(Δ, 0)] .= -1/3 * π^2 / h^2\n  for i = 2:n\n    Δ[diagind(Δ, i - 1)] = \n    Δ[diagind(Δ, 1 - i)] .= 2 * (-1)^i / (i - 1)^2 / h^2\n  end\n  return Δ\nend\n\n# make logarithmic grid on r\nx = -30.0:0.175:5.0; r = exp.(x); r² = exp.(2x)\n\n# set up parameters of the hydrogen atom hamiltonian:\n# reduced mass, azimuthal quantum number, and nucleus charge\nμ = 1; l = 0; Z = 1;\n\n# make hamiltonian matrix\nH = -1/2μ*laplacian(x) + Diagonal(1/2μ * (l + 1/2)^2 .- Z*r);\n\n# form positive definite well-conditioned matrix S\nα = 1e5; β = 1;\nS = β*H + α * Diagonal(r .* r);\n\n# solve generalized eigenproblem\nθ, y = eigen(H, S);\nϵ = α * θ ./ (1.0 .- θ);\n\n# see energy of the 1s-state of hydrogen atom\nϵ[1]","category":"page"},{"location":"#Functions-for-solving-Schrödinger-equation","page":"AtomEnergyLevels.jl","title":"Functions for solving Schrödinger equation","text":"","category":"section"},{"location":"#Differentiation-matrix","page":"AtomEnergyLevels.jl","title":"Differentiation matrix","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"laplacian","category":"page"},{"location":"#AtomEnergyLevels.laplacian","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.laplacian","text":"laplacian(n::Int, h::Real)\nlaplacian(x::AbstractRange{T}) where {T}\n\nReturns n times n differentiation matrix mathrmD^(2) using sinc interpolants.  See Weideman, J.A., and Reddy, S.C. (2000). A MATLAB differentiation matrix suite.  ACM Transactions on Mathematical Software (TOMS), 26(4), 465-519. Equation (20) on page 485.\n\nOn input: n - size of a uniform grid with step size h, or x - AbstractRange object.\n\nExample\n\nSecond derivative of the gaussian function\n\n    partial^2 exp(-x^2)partial x^2 = 2(2x^2 - 1)exp(-x^2)\n\njulia> x = -6:0.01:6\n-6.0:0.01:6.0\n\njulia> d²f_exact = 2(2x.^2 .- 1).*exp.(-x.^2);\n\njulia> d²f_approx = laplacian(length(x), step(x)) * exp.(-x.^2);\n\njulia> isapprox(d²f_exact, d²f_approx, atol = 1e-10)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Using the following code, one can find the vibrational levels of the radical OH⋅, for which the potential energy surface is approximated through the Morse potential.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"using AtomEnergyLevels\nusing LinearAlgebra, Test, Printf\n\n# Morse potential\nV(x, D, β, x₀) = D*(exp(-β*(x-x₀))-1)^2 - D;\n\n# parameters of the O-H bond\nD  = 0.1994;   β = 1.189;  x₀ = 1.821;\nmH = 1.00794; mO = 15.9994; μ = 1822.8885*(mH*mO)/(mH+mO);  \n\n# grid\nx = (2.0 - x₀):0.1:(12.0 + x₀);\n\n# hamiltonian\nH = -1/2μ*laplacian(x) + Diagonal(V.(x, D, β, x₀));\n\n# numerical solution\nϵ, ψ = eigen(H);                      \n\n# exact solution\nω = β*sqrt(2D/μ); δ = ω^2 / 4D;  \nE(n) = ω*(n+1/2) - δ*(n+1/2)^2 - D;\n\n# comparison\nfor i in 1:5\n  @printf \"Level %i: E(exact) = %5.10f E(approx) = %5.10f\\n\" i E(i-1) ϵ[i]\nend","category":"page"},{"location":"#Atomic-electron-configuration-parsing","page":"AtomEnergyLevels.jl","title":"Atomic electron configuration parsing","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"conf_enc","category":"page"},{"location":"#AtomEnergyLevels.conf_enc","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.conf_enc","text":"conf_enc(input::AbstractString; maxn = 7, spins = 2,\nnotation = Dict(\"s\" => 0, \"p\" => 1, \"d\" => 2, \"f\" => 3))\n\nParse an input string as the configuration of an atom. An electronic  configuration is string like \"[He] 2s2 2p1\", where the notation [X] indicates  that all subshells associated with the noble gas X are fully occupied.  The function encodes the configuration provided \"[He] 2s2 2p1\" and returns  ((2.0, 2.0), (1.0,)), where subtuples (2.0, 2.0) and (1.0,) are  populations of all the l = 0 and l = 1 shells, i.e. ((1s², 2s²), (2p¹,)).  Within each subtuple populations of levels with the same azimuthal quantum  number l and subsequent radial quantum numbers nᵣ are listed.\n\nOptional keyword arguments are:\n\nmaxn  - maximal permitted principal quantum number\nspins - number of spin up/spin down electrons per level (1 or 2)\nnotation - dictionary encoding the correspondence of quantum numbers l = 0, 1... to letters: s, p, d... See spectroscopic notation\n\nExamples\n\njulia> conf = conf_enc(\"[Ar] 3d1 4s2\") # Scandium, Z = 21\n((2.0, 2.0, 2.0, 2.0), (6.0, 6.0), (1.0,))\n\njulia> sum(collect(Iterators.flatten(conf)))\n21.0\n\n\n\n\n\n","category":"function"},{"location":"#Radial-Schrödinger-equation-solver","page":"AtomEnergyLevels.jl","title":"Radial Schrödinger equation solver","text":"","category":"section"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"radial_shr_eq","category":"page"},{"location":"#AtomEnergyLevels.radial_shr_eq","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.radial_shr_eq","text":"radial_shr_eq(V::Function = r -> -1/r, x::AbstractRange{T} = -30.0:0.1:5.0;\nconf = 1, μ = 1.0, α = 1e5) where {T}\n\nradial_shr_eq(V::Array{S}, x::AbstractRange{T}; \nconf = 1, μ = 1.0, α = 1e5) where {S, T}\n\nSolve Schrödinger equation (atomic units are assumed)\n\n-frac12muleftfrac1r^2fracpartialpartial rleft(r^2fracpartial R_nlpartial rright)-fracl(l+1)r^2R_nl(r)right+V(r)R_nl(r)=E_nlR_nl(r)\n\nfor the particle in a spherically symmetric potential.\n\nOn input:\n\nV - potential, provided as an explicit function V(r) (by default r -> -1/r), or an array of values on a radial grid.\nx - uniform grid (AbstractRange object), such as r_i = exp(x_i). Default grid is -30.0:0.1:5.0\nconf - electronic configuration of interest, by default: ground state conf = ((1.0),)\nμ - reduced mass, by default μ = 1.0 (electron mass in a.u.)\nα - parameter α = 10^5 is used in the matrix pencil mathrmH + αmathrmS for the generalized eigenproblem.\n\nOn output function returns:\n\nsum of one-particle energies: E = sum_i=1^N n_i epsilon_i\nparticles density: rho(x) = frac14pi sum_i=1^N n_i y_i^2(x)\norbitals y_i(x), corresponding eigenvalues epsilon_i (energy levels), azimuthal and radial quantum numbers l, n_r, and level populations n_i (as listed in conf).\n\nExamples\n\nHydrogen atom\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"1s1\")).energy ≈ -1/2\ntrue\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"2s1\")).energy ≈ -1/8\ntrue\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"2p1\")).energy ≈ -1/8\ntrue\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"3s1\")).energy ≈ -1/18\ntrue\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"3p1\")).energy ≈ -1/18\ntrue\n\njulia> radial_shr_eq(r -> -1/r, conf = conf_enc(\"3d1\")).energy ≈ -1/18\ntrue\n\n\n\n\n\n","category":"function"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"Using the following code, one can find the energy levels for a spherically-symmetric three-dimensional harmonic oscillator.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"To find the energy levels with quantum numbers nᵣ = 0, 1, 2 and l = 0, 1, 2, the levels of interest should be included in the electronic configuration.","category":"page"},{"location":"","page":"AtomEnergyLevels.jl","title":"AtomEnergyLevels.jl","text":"using AtomEnergyLevels, Printf\n\nfunction isotropic_harmonic_oscillator(cfg)\n  ψ = radial_shr_eq(r -> 1/2*r^2, conf = conf_enc(cfg)).orbitals;\n  @printf(\"\\tϵ(calc.)\\tϵ(exact)\\tΔϵ\\n\");\n  for (quantum_numbers, orbital) in sort(collect(ψ), by = x -> last(x).ϵᵢ)\n    nᵣ, l = quantum_numbers\n    ϵ_calc = orbital.ϵᵢ\n    n = nᵣ + l + 1\n    ϵ_exact = 2nᵣ + l + 3/2\n    @printf(\"%i%s\\t%10.8f\\t%10.8f\\t%+0.6e\\n\",\n            n, atomic_shell[l], ϵ_calc, ϵ_exact, ϵ_exact - ϵ_calc)\n  end\nend\n\nisotropic_harmonic_oscillator(\"1s1 2s1 3s1 2p1 3p1 4p1 3d1 4d1 5d1\");","category":"page"}]
}
