"""
    laplacian(n::Int, h::Real)
    laplacian(x::AbstractRange{T}) where {T}

Returns ``n \\times n`` differentiation matrix ``\\mathrm{D^{(2)}}`` using sinc interpolants. 
See Weideman, J.A., and Reddy, S.C. (2000). A MATLAB differentiation matrix suite. 
ACM Transactions on Mathematical Software (TOMS), 26(4), 465-519.
Equation (20) on page 485.

On input: `n` - size of a uniform grid with step size `h`, or `x` - `AbstractRange` object.

# Example
Second derivative of the gaussian function
```math 
    \\partial^2 \\exp(-x^2)/\\partial x^2 = 2(2x^2 - 1)\\exp(-x^2)
```
```jldoctest
julia> x = -6:0.01:6
-6.0:0.01:6.0

julia> d²f_exact = 2(2x.^2 .- 1).*exp.(-x.^2);

julia> d²f_approx = laplacian(length(x), step(x)) * exp.(-x.^2);

julia> isapprox(d²f_exact, d²f_approx, atol = 1e-10)
true
```
"""
laplacian(x::AbstractRange{T}) where {T} = _laplacian(length(x), step(x))
laplacian(n::Int, h::Real) = _laplacian(n, h)

function _laplacian(n, h)
  Δ = Matrix{Float64}(undef, n, n)
  Δ[diagind(Δ, 0)] .= -1/3*π^2 / h^2
  for i=2:n
    @inbounds Δ[diagind(Δ, i-1)] = Δ[diagind(Δ, 1-i)] .= 2*(-1)^i / (i-1)^2 / h^2
  end
  return Δ
end
