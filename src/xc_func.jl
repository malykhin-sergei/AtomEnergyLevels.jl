@enum XC_Type begin
  xc_type_lda
  xc_type_gga
  xc_type_mgga
  xc_type_hybrid
end

struct XC_Functional 
  kind::XC_Type
  xc_func::Function
end

function xc_func(kind::XC_Type = xc_type_lda, spin_polarization::Bool = false, 
                 exchange::Function = LDA_X, correlation::Function = LDA_C_VWN, threshold = eps())
  function xc_lda(ρ, vxc = similar(ρ), exc = similar(ρ))
    if !(length(ρ) == length(vxc) == length(exc))
      throw(DimensionMismatch("Dimensions of the arrays ρ, vxc, and exc doesn't match"))
    end
    @inbounds for i in eachindex(ρ)
      vxc[i], exc[i] = ρ[i] ≤ threshold ? (0, 0) : exchange(ρ[i]) .+ correlation(ρ[i])
    end
    return vxc, exc
  end
  
  if kind == xc_type_lda && spin_polarization == false
    return XC_Functional(kind, xc_lda)
  else
    return nothing
  end
end
