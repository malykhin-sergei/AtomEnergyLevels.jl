module AtomEnergyLevels

import LinearAlgebra: Symmetric, eigen, diagind, I 
import Printf: @sprintf

include("integrate.jl")
include("sincdif.jl")
include("conf_parse.jl")
include("periodic_table.jl")
include("xclib.jl")
include("xc_func.jl")
include("rshreq.jl")
include("TF.jl")
include("lda.jl")

export laplacian, radial_shr_eq, TF, lda
export LDA_X, LDA_C_CHACHIYO, LDA_C_VWN
export XC_Type, XC_Functional, xc_func
export xc_type_lda, xc_type_gga, xc_type_mgga, xc_type_hybrid
export conf_enc, Nₑ, @c_str
export shells, atom
export laplacian
export radial_shr_eq

end # module
