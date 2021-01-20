module AtomEnergyLevels

include("integrate.jl")
include("sincdif.jl")
include("conf_parse.jl")
include("periodic_table.jl")
include("dft_xc_functionals.jl")
include("rshreq.jl")
include("TF.jl")
include("lda.jl")

import LinearAlgebra: Symmetric, eigen, diagind, I 
import Printf: @sprintf

export laplacian, radial_shr_eq, TF, lda
export LDA_X, LDA_C_CHACHIYO, LDA_C_VWN
export SVWN!, Xα!
export conf_enc, Nₑ, @c_str
export shells, atom
export laplacian
export radial_shr_eq

end # module
