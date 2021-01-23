"""
    LDA_X(ρ; α = 2/3)

Gives Slater's approximation of the exchange energy density 
``\\epsilon_x``, and potential ``V_x``.
    
## References

1. P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
2. F. Bloch, Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
3. J.C. Slater, Phys. Rev. 81, 385 (1951) (doi: 10.1103/PhysRev.81.385)
"""
function LDA_X(ρ; α = 2/3)
  vx = -3α * (3ρ/8π)^(1/3)
  ϵx = 3/4 * vx
  return vx, ϵx
end

"""
    LDA_C_XALPHA(ρ; α = 1)
    
Slater’s Xα

## References

1. J.C. Slater, Phys. Rev. 81, 385 (1951) (doi: 10.1103/PhysRev.81.385)
2. J.C. Slater, T.M. Wilson, and J.H. Wood, Phys. Rev. 179, 28 (1969)
"""
LDA_C_XALPHA(ρ; α = 1) = LDA_X(ρ, α = α - 2/3)

"""
    LDA_C_CHACHIYO(ρ; a  = -0.01554535, b₁ = 20.4562557, b  = 20.4562557)

Gives Chachiyo's approximation of the correlation energy density 
``\\epsilon_x``, and potential ``V_x``.
    
## References

1. T. Chachiyo, J. Chem. Phys. 145, 021101 (2016) (doi: 10.1063/1.4958669)
"""
function LDA_C_CHACHIYO(ρ; a  = -0.01554535, b₁ = 20.4562557, b  = 20.4562557)
  rs = (3/4π / ρ)^(1/3)
  ϵc = a * log(1 + b₁ / rs + b / rs^2)
  vc = ϵc + a * (b₁ * rs + 2b) / (3 * (rs^2 + b₁ * rs + b))
  return vc, ϵc
end

"""
    LDA_C_VWN(ρ; a = 0.0310907, b = 3.72744, c = 12.9352, x0 = -0.10498)

Gives VWN5 approximation of the correlation energy density 
``\\epsilon_x``, and potential ``V_x``.
    
## References

1. S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
"""
function LDA_C_VWN(ρ; a = 0.0310907, b = 3.72744, c = 12.9352, x0 = -0.10498)
  rs = (3/4π/ρ)^(1/3)
  q = sqrt(4.0 * c - b * b)
  f1 = 2.0 * b / q
  f2 = b * x0 / (x0 * x0 + b * x0 + c)
  f3 = 2.0 * (2.0 * x0 + b) / q
  rs12 = sqrt(rs)
  fx = rs + b * rs12 + c
  qx = atan(q / (2.0 * rs12 + b) )
  ϵc = a * (log(rs/fx) + f1 * qx - f2 * (log((rs12 - x0)^2 / fx) + f3 * qx))
  tx = 2.0 * rs12 + b
  tt = tx * tx + q * q
  vc = ϵc - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b / tt - f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b) / tt))
  return vc, ϵc
end

function xc_lda(exch_fun, corr_fun, threshold = eps())
  function xc!(ρ, vxc, exc)
    @simd for i in eachindex(ρ)
      @inbounds vxc[i], exc[i] = ρ[i] ≤ threshold ? (0, 0) : exch_fun(ρ[i]) .+ corr_fun(ρ[i])
    end
    return vxc, exc
  end
end

# TODO:
# V. V. Karasiev, J. Chem. Phys. 145, 157101 (2016) (doi: 10.1063/1.4964758)
# LDA_C_KARASIEV_MOD(ρ) = LDA_C_CHACHIYO(ρ, b₁ = 21.7392245)