
struct Data_CUDA
    ID_sp_dest:: Array
    Married:: Array
    MarriedHH:: CuArray
    OriginalSpouse_:: CuArray
    Sex_:: CuArray
    insample:: CuArray
    insample!_Bool:: CuArray
    firstsampled:: CuArray
    firstsampled_Bool:: CuArray
    I_st1:: CuArray
    I_st2_single:: CuArray
    I_st2_cpl:: CuArray
    I_st1_:: CuArray
    I_st2_single_:: CuArray
    I_st2_cpl_:: CuArray
    marriedhh:: CuArray
    married_tminus1_:: CuArray
    sex_:: CuArray
    married_:: CuArray
    originalspouse:: CuArray
    originalspouse_reduced:: CuArray
    First_t_with_sp:: CuArray
    same_sp:: CuArray
    NotMarried:: CuArray
    NotMarried_sort:: CuArray
    OriginalSpouse_sort:: CuArray
    insample_sort:: CuArray
    idx_NotMarried_sort:: CuArray
    idx_OriginalSpouse_sort:: CuArray
    idx_insample_sort:: CuArray
    idx_p_st1_m_singleb4:: CuArray
    idx_p_st1_f_singleb4:: CuArray
    idx_p0_st1_m_singleb4:: CuArray
    idx_p0_st1_f_singleb4:: CuArray
    idx_p_st1_m_marriedb4:: CuArray
    idx_p_st1_f_marriedb4:: CuArray
    idx_p0_st1_m_marriedb4:: CuArray
    idx_p0_st1_f_marriedb4:: CuArray
    idx_p_st2_m_single:: CuArray
    idx_p_st2_f_single:: CuArray
    idx_p_st2_cpl:: CuArray
    idx_p_st1_m_singleb4_reduced:: CuArray
    idx_p_st1_f_singleb4_reduced:: CuArray
    idx_p0_st1_m_singleb4_reduced:: CuArray
    idx_p0_st1_f_singleb4_reduced:: CuArray
    idx_p_st1_m_marriedb4_reduced:: CuArray
    idx_p_st1_f_marriedb4_reduced:: CuArray
    idx_p0_st1_m_marriedb4_reduced:: CuArray
    idx_p0_st1_f_marriedb4_reduced:: CuArray
    idx_p_st2_m_single_reduced:: CuArray
    idx_p_st2_f_single_reduced:: CuArray
    idx_p_st2_cpl_reduced:: CuArray
    partidx_p_st1_m_singleb4_data:: CuArray
    partidx_p_st1_m_marriedb4_data:: CuArray
    partidx_p_st1_f_singleb4_data:: CuArray
    partidx_p_st1_f_marriedb4_data:: CuArray
    partidx_p_st2_f_data:: CuArray
    partidx_p_st2_m_data:: CuArray
    partidx_p_st2_cpl_data:: CuArray
    idx_p_st1_m_singleb4_reduced_obs:: CuArray
    idx_p_st1_f_singleb4_reduced_obs:: CuArray
    idx_p0_st1_m_singleb4_reduced_obs:: CuArray
    idx_p0_st1_f_singleb4_reduced_obs:: CuArray
    idx_p_st1_m_marriedb4_reduced_obs:: CuArray
    idx_p_st1_f_marriedb4_reduced_obs:: CuArray
    idx_p0_st1_m_marriedb4_reduced_obs:: CuArray
    idx_p0_st1_f_marriedb4_reduced_obs:: CuArray
    idx_p_st2_m_single_reduced_obs:: CuArray
    idx_p_st2_f_single_reduced_obs:: CuArray
    idx_p_st2_cpl_reduced_obs:: CuArray

    idx_p_st1_m_singleb4_choicesonly:: CuArray
    idx_p_st1_m_marriedb4_choicesonly:: CuArray


    t:: CuArray
    phistate_sp_orig_:: CuArray
    hstate_sp_orig_:: CuArray
    s_sp_orig_:: CuArray
    phistate_:: CuArray
    hstate_:: CuArray
    s_:: CuArray
    xstate_:: CuArray
    xstate_m_cpl_:: CuArray
    hstate_m_cpl_:: CuArray
    s_m_cpl_:: CuArray
    xstate_f_cpl_:: CuArray
    hstate_f_cpl_:: CuArray
    s_f_cpl_:: CuArray
end


struct Data_redform_CUDA
    xx_st1:: CuArray
    xx_st2:: CuArray
    xx_st2_cpl:: CuArray
    dim_x_st1_2:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_x_st2_2:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_x_st2_cpl_2:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_y_st1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_y_st2:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_y_st2_cpl:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_x_st1_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_x_st2_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_x_st2_cpl_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    originalspouse:: CuArray
    originalspouse_reduced:: CuArray
    
    idx_st1_m_allchoices:: CuArray
    idx_st1_m_singleb4_allchoices_observedstate:: CuArray
    idx_st1_m_marriedb4_allchoices_observedstate:: CuArray
end

struct Data_grid_CUDA
    RX1_st1_m_marriedb4:: CuArray
    RX1_st1_f_marriedb4:: CuArray
    RX1_st1_m_singleb4:: CuArray
    RX1_st1_f_singleb4:: CuArray
    RX1_st2_m:: CuArray
    RX1_st2_f:: CuArray
    RX1_st2_cpl:: CuArray
end

struct XstructGrid_CUDA
    Xstruct_st1_marriedb4:: CuArray
    Xstruct_st1_singleb4:: CuArray
    Xstruct_st2_m_married:: CuArray
    Xstruct_st2_f_married:: CuArray
    Xstruct_st2_single:: CuArray
    idx_st2_into_st1_f:: CuArray
    idx_st2_into_st1_m:: CuArray

    dim_Xstruct_st1_marriedb4_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_Xstruct_st1_singleb4_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_Xstruct_st2_single_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_Xstruct_st2_m_married_1:: Union{BigInt, Int64, Int32, Int16, Int8}
    dim_Xstruct_st2_f_married_1:: Union{BigInt, Int64, Int32, Int16, Int8}
end

struct CCP_True_CUDA
    p0_st1_m_marriedb4:: CuArray
    p0_st1_f_marriedb4:: CuArray
    p0_st1_m_singleb4:: CuArray
    p0_st1_f_singleb4:: CuArray
    p0_st2_m:: CuArray
    p0_st2_f:: CuArray
    p0_st2_cpl:: CuArray
    p_st1_m_marriedb4:: CuArray
    p_st1_f_marriedb4:: CuArray
    p_st1_m_singleb4:: CuArray
    p_st1_f_singleb4:: CuArray
    p_st2_m:: CuArray
    p_st2_f:: CuArray
    p_st2_cpl:: CuArray
end

struct Data_Struct_GMM_CUDA
    Z_st1_m_singleb4_data:: CuArray
    Z_st1_m_marriedb4_data:: CuArray
    Z_st1_f_singleb4_data:: CuArray
    Z_st1_f_marriedb4_data:: CuArray
    Z_st2_m_data:: CuArray
    Z_st2_f_data:: CuArray
    Z_st2_cpl_data:: CuArray
end

function adapt_structure_to_gpu(dt_CPU,dt_rf_CPU,dt_rf_grid_CPU, Xstruct_grid_CPU, CCP_true_CPU, dt_z_GMM_CPU)
    #CUDA doesn't handle missing values, convert to NaN
    # also, coalesce converts data type to Real, which is not compatible with CUDA, so I have to convert it back
    dt_CUDA = Data_CUDA(
    dt_CPU.ID_sp_dest,
    dt_CPU.Married,
    CUDA.adapt(CuArray, dt_CPU.MarriedHH),
    CUDA.adapt(CuArray, dt_CPU.OriginalSpouse_),
    CUDA.adapt(CuArray, dt_CPU.Sex_),
    CUDA.adapt(CuArray, dt_CPU.insample),
    CUDA.adapt(CuArray, dt_CPU.insample!_Bool),
    CUDA.adapt(CuArray, dt_CPU.firstsampled),
    CUDA.adapt(CuArray, dt_CPU.firstsampled_Bool),
    CUDA.adapt(CuArray, dt_CPU.I_st1),
    CUDA.adapt(CuArray, dt_CPU.I_st2_single),
    CUDA.adapt(CuArray, dt_CPU.I_st2_cpl),
    CUDA.adapt(CuArray, dt_CPU.I_st1_),
    CUDA.adapt(CuArray, dt_CPU.I_st2_single_),
    CUDA.adapt(CuArray, dt_CPU.I_st2_cpl_),
    CUDA.adapt(CuArray, dt_CPU.marriedhh),
    CUDA.adapt(CuArray, dt_CPU.married_tminus1_),
    CUDA.adapt(CuArray, dt_CPU.sex_),
    CUDA.adapt(CuArray, dt_CPU.married_),
    CUDA.adapt(CuArray, dt_CPU.originalspouse),
    CUDA.adapt(CuArray, dt_CPU.originalspouse_reduced),
    CUDA.adapt(CuArray, dt_CPU.First_t_with_sp),
    CUDA.adapt(CuArray, dt_CPU.same_sp),
    CUDA.adapt(CuArray, dt_CPU.NotMarried),
    CUDA.adapt(CuArray, dt_CPU.NotMarried_sort),
    CUDA.adapt(CuArray, dt_CPU.OriginalSpouse_sort),
    CUDA.adapt(CuArray, dt_CPU.insample_sort),
    CUDA.adapt(CuArray, dt_CPU.idx_NotMarried_sort),
    CUDA.adapt(CuArray, dt_CPU.idx_OriginalSpouse_sort),
    CUDA.adapt(CuArray, dt_CPU.idx_insample_sort),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_singleb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_singleb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_singleb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_singleb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_marriedb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_marriedb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_marriedb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_marriedb4),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_m_single),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_f_single),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_cpl),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_singleb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_singleb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_singleb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_singleb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_marriedb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_marriedb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_marriedb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_marriedb4_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_m_single_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_f_single_reduced),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_cpl_reduced),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st1_m_singleb4_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st1_m_marriedb4_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st1_f_singleb4_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st1_f_marriedb4_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st2_f_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st2_m_data),
    CUDA.adapt(CuArray, dt_CPU.partidx_p_st2_cpl_data),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_singleb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_singleb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_singleb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_singleb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_marriedb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_f_marriedb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_m_marriedb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p0_st1_f_marriedb4_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_m_single_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_f_single_reduced_obs),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st2_cpl_reduced_obs),
    

    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_singleb4_choicesonly),
    CUDA.adapt(CuArray, dt_CPU.idx_p_st1_m_marriedb4_choicesonly),
    
    CUDA.adapt(CuArray, dt_CPU.t),
    CUDA.adapt(CuArray, dt_CPU.phistate_sp_orig_),
    CUDA.adapt(CuArray, dt_CPU.hstate_sp_orig_),
    CUDA.adapt(CuArray, dt_CPU.s_sp_orig_),
    CUDA.adapt(CuArray, dt_CPU.phistate_),
    CUDA.adapt(CuArray, dt_CPU.hstate_),
    CUDA.adapt(CuArray, dt_CPU.s_),
    CUDA.adapt(CuArray, dt_CPU.xstate_),
    CUDA.adapt(CuArray, dt_CPU.xstate_m_cpl_),
    CUDA.adapt(CuArray, dt_CPU.hstate_m_cpl_),
    CUDA.adapt(CuArray, dt_CPU.s_m_cpl_),
    CUDA.adapt(CuArray, dt_CPU.xstate_f_cpl_),
    CUDA.adapt(CuArray, dt_CPU.hstate_f_cpl_),
    CUDA.adapt(CuArray, dt_CPU.s_f_cpl_))

    #CUDA doesn't handle missing values, convert to NaN
    # also, coalesce converts data type to Real, which is not compatible with CUDA, so I have to convert it back
    dt_rf_CUDA=Data_redform_CUDA(
    CUDA.adapt(CuArray, precision.(coalesce.(dt_rf_CPU.xx_st1,0))),
    CUDA.adapt(CuArray, precision.(coalesce.(dt_rf_CPU.xx_st2,0))),
    CUDA.adapt(CuArray, precision.(coalesce.(dt_rf_CPU.xx_st2_cpl,0))),
    dt_rf_CPU.dim_x_st1_2,
    dt_rf_CPU.dim_x_st2_2,
    dt_rf_CPU.dim_x_st2_cpl_2,
    dt_rf_CPU.dim_y_st1,
    dt_rf_CPU.dim_y_st2,
    dt_rf_CPU.dim_y_st2_cpl,
    dt_rf_CPU.dim_x_st1_1,
    dt_rf_CPU.dim_x_st2_1,
    dt_rf_CPU.dim_x_st2_cpl_1,
    CUDA.adapt(CuArray, dt_rf_CPU.originalspouse_reduced),
    CUDA.adapt(CuArray, dt_rf_CPU.originalspouse),

    CUDA.adapt(CuArray, dt_rf_CPU.idx_st1_m_allchoices),
    CUDA.adapt(CuArray, dt_rf_CPU.idx_st1_m_singleb4_allchoices_observedstate),
    CUDA.adapt(CuArray, dt_rf_CPU.idx_st1_m_marriedb4_allchoices_observedstate))

    dt_rf_grid_CUDA= Data_grid_CUDA(
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st1_m_marriedb4), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st1_f_marriedb4), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st1_m_singleb4), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st1_f_singleb4), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st2_m), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st2_f), 
    CUDA.adapt(CuArray, dt_rf_grid_CPU.RX1_st2_cpl));

    Xstruct_grid_CUDA=XstructGrid_CUDA(
        CUDA.adapt(CuArray, Xstruct_grid_CPU.Xstruct_st1_marriedb4),
        CUDA.adapt(CuArray, Xstruct_grid_CPU.Xstruct_st1_singleb4),
        CUDA.adapt(CuArray, Xstruct_grid_CPU.Xstruct_st2_m_married),
        CUDA.adapt(CuArray, Xstruct_grid_CPU.Xstruct_st2_f_married),
        CUDA.adapt(CuArray, Xstruct_grid_CPU.Xstruct_st2_single) ,
        CUDA.adapt(CuArray, Xstruct_grid_CPU.idx_st2_into_st1_f) ,
        CUDA.adapt(CuArray, Xstruct_grid_CPU.idx_st2_into_st1_m), 

        Xstruct_grid_CPU.dim_Xstruct_st1_marriedb4_1 ,
        Xstruct_grid_CPU.dim_Xstruct_st1_singleb4_1,
        Xstruct_grid_CPU.dim_Xstruct_st2_single_1 ,
        Xstruct_grid_CPU.dim_Xstruct_st2_m_married_1,
        Xstruct_grid_CPU.dim_Xstruct_st2_f_married_1
    )

    CCP_true_CUDA=CCP_True_CUDA(
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st1_m_marriedb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st1_f_marriedb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st1_m_singleb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st1_f_singleb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st2_m),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st2_f),
        CUDA.adapt(CuArray, CCP_true_CPU.p0_st2_cpl),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st1_m_marriedb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st1_f_marriedb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st1_m_singleb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st1_f_singleb4),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st2_m),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st2_f),
        CUDA.adapt(CuArray, CCP_true_CPU.p_st2_cpl)
    )

    data_struct_GMM_CUDA=Data_Struct_GMM_CUDA(
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st1_m_singleb4_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st1_m_marriedb4_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st1_f_singleb4_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st1_f_marriedb4_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st2_m_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st2_f_data),
        CUDA.adapt(CuArray, dt_z_GMM_CPU.Z_st2_cpl_data)
    )



    return dt_CUDA, dt_rf_CUDA, dt_rf_grid_CUDA, Xstruct_grid_CUDA, CCP_true_CUDA, data_struct_GMM_CUDA
end

