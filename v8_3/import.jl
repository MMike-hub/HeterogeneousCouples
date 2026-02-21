include("data_genbus4_solve.jl")
include("data_genbus4_simulate.jl")
include("data_generate_wrap.jl")
include("data_compress.jl")
include("data_adapt_structure_to_gpu.jl")
include("data_struct_stategrid.jl")
include("data_struct_GMM.jl")
include("placeholders.jl")

include("ForwardValues.jl")

include("CCP_structural.jl")
include("CCP_reducedform.jl")
include("likeCCP_rf.jl")
include("likeCCP_struct.jl")

include("reducedform_wrap.jl")
include("reducedform_gmm.jl")
include("reducedform_pseudoMLE.jl")

include("struct_wrap.jl")
include("struct_pseudoMLE.jl")
include("struct_GMM.jl")

include("xgrid.jl")
include("gmm_struct.jl")
include("fhat_gmm.jl")
include("data_rf_state_original.jl")
include("data_rf_state_wrap.jl")

include("resCCP.jl")
include("get_slice.jl")

include("update_ptype_kappanoorigsp.jl")
include("update_ptype_kappanoorigsp_old.jl")
include("paramsetup.jl")