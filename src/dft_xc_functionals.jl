function Xα(ρ; α = 2/3)
  # J.C. Slater, Phys. Rev. 81, 385 (1951)
  vₓ = -3α * (3ρ/8π)^(1/3)
  ϵₓ = 3/4 * vₓ
  return vₓ, ϵₓ
end

function SVWN(ρ)
  ρ <= eps() && return 0.0, 0.0
  vₓ, ϵₓ = Xα(ρ)
  # S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
  rs = (3/4π/ρ)^(1/3)
  a = 0.0310907; b = 3.72744; c = 12.9352; x0 = -0.10498;
  q = sqrt(4.0 * c - b * b)
  f1 = 2.0 * b / q
  f2 = b * x0 / (x0 * x0 + b * x0 + c)
  f3 = 2.0 * (2.0 * x0 + b) / q
  rs12 = sqrt(rs)
  fx = rs + b * rs12 + c
  qx = atan(q / (2.0 * rs12 + b) )
  ec = a * (log(rs/fx) + f1 * qx - f2 * (log((rs12 - x0)^2 / fx) + f3 * qx))
  tx = 2.0 * rs12 + b
  tt = tx * tx + q * q
  vc = ec - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b /
     tt - f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b) / tt))
  return vₓ + vc, ϵₓ + ec
end
