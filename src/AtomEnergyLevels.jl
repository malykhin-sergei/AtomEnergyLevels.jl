module AtomEnergyLevels

import LinearAlgebra: Symmetric, eigen, diagind, I 
import Printf: @sprintf

include("integrate.jl")
include("sincdif.jl")
include("conf_parse.jl")
include("periodic_table.jl")
include("dft_xc_functionals.jl")
include("rshreq.jl")
include("TF.jl")
include("lda.jl")

export laplacian, radial_shr_eq, TF, lda
export xc_lda, LDA_X, LDA_C_XALPHA, LDA_C_CHACHIYO, LDA_C_VWN
export conf_enc, Nₑ, @c_str
export shells, atom
export laplacian
export radial_shr_eq

end # module
