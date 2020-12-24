"""
    TF(r, Z)

Approximate solution of the [Thomas-Fermi equation](https://en.wikipedia.org/wiki/Thomas-Fermi_model) 
for neutral atom with nuclear charge `Z`. Returns density at the distance `r`.

See Moliere, G. (1947). Theorie der streuung schneller geladener teilchen i. einzelstreuung am abgeschirmten coulomb-feld.
Zeitschrift für Naturforschung A, 2(3), 133-145.

## Example
For a density function ``\\rho(r)`` of ``N`` particles the following equation holds.
```math
N = 4\\pi\\int dr \\rho(r)r^2
``` 
On the grid ``r_i = \\exp(x_i)`` it is
```math 
N = 4\\pi\\int dx \\rho(r(x))r^3
```
```jldoctest
julia> x = -30:0.1:20
-30.0:0.1:20.0

julia> r, n, dx = exp.(x), length(x), step(x);

julia> isapprox(4π * dx * sum(TF.(r, 18) .* r .^3), 18, atol=1e-10)
true
```
"""
function TF(r, Z)
  b = 1/2*(3π/4)^(2/3)*Z^(-1/3)
  x = r/b

  B₁, B₂, B₃ = 0.1, 0.55, 0.35
  β₁, β₂, β₃ = 6.0, 1.20, 0.30

  Z/(4π*b^3*x)*(B₁*β₁^2*exp(-β₁*x) + B₂*β₂^2*exp(-β₂*x) + B₃*β₃^2*exp(-β₃*x))
end
