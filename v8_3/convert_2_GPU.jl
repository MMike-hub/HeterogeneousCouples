function f_convert_2_GPU()
    global dt,  Xstruct_grid, CCP_true, dt_struct_GMM, PType, PType_true,
    xtran, xtran_married, bccp, bccp_true, Lambda, Pi_cpl, Pi_single, Pi_single_m, Pi_single_f,  xval, hval, 
    zval, alpha_wage, bwage, 
    dt_CPU, dt_rf_CPU, dt_rf_grid_CPU, Xstruct_grid_CPU, CCP_true_CPU, dt_struct_GMM_CPU, xval_CPU, zval_CPU, hval_CPU,
    b1_CPU, PType_CPU, PType_true_CPU, xtran_CPU, xtran_married_CPU, bccp_CPU, bccp_true_CPU, Lambda_CPU, Pi_cpl_CPU, Pi_single_m_CPU, Pi_single_f_CPU, alpha_wage_CPU, bwage_CPU, Pi_single_CPU

    global zeros_dim_x_st2_2, dummy_size0, dummy_size10, dummy_sizeNT, 
    zeros_st1_m_part1, zeros_st1_m_part2, zeros_st1_f_part1, zeros_st1_f_part2, zeros_st2_cpl_m, zeros_st2_cpl_f, zeros_st2_m, 
    zeros_st2_f, zeros_dim_x_st2_cpl_2, ones_st2_married, ones_st2_single
    # Convert arrays in structures to CuArray
        dt_CPU=dt
        dt_rf_CPU=dt_rf
        dt_rf_grid_CPU=dt_rf_grid
        Xstruct_grid_CPU=Xstruct_grid
        CCP_true_CPU=CCP_true
        dt_struct_GMM_CPU=dt_struct_GMM

        
        dt_CUDA, dt_rf_CUDA, dt_rf_grid_CUDA, Xstruct_grid_CUDA, CCP_true_CUDA, 
        dt_struct_GMM_CUDA=adapt_structure_to_gpu(dt_CPU, dt_rf_CPU, dt_rf_grid_CPU, Xstruct_grid_CPU, CCP_true_CPU, dt_struct_GMM_CPU);
        
        dt=dt_CUDA
        dt_rf=dt_rf_CUDA
        dt_rf_grid=dt_rf_grid_CUDA
        Xstruct_grid=Xstruct_grid_CUDA
        CCP_true=CCP_true_CUDA
        dt_struct_GMM=dt_struct_GMM_CUDA


        # Convert parameters vectors to CuArray
        b1_CPU=b1
        PType_CPU=PType
        PType_true_CPU=PType_true
        xtran_CPU=xtran
        xtran_married_CPU=xtran_married
        bccp_CPU=bccp
        bccp_true_CPU=bccp_true
        alpha_wage_CPU=alpha_wage
        bwage_CPU=bwage
        Lambda_CPU=Lambda
        Pi_cpl_CPU=Pi_cpl
        Pi_single_m_CPU=Pi_single_m
        Pi_single_f_CPU=Pi_single_f
        Pi_single_CPU=Pi_single
        xval_CPU=xval
        zval_CPU=zval
        hval_CPU=hval

        b1=CUDA.adapt(CuArray,b1)
        PType=CUDA.adapt(CuArray,PType)
        PType_true=CUDA.adapt(CuArray,PType_true)
        xtran=CUDA.adapt(CuArray, xtran)
        xtran_married=CUDA.adapt(CuArray, xtran_married)
        bccp=CUDA.adapt(CuArray, bccp)
        bccp_true=CUDA.adapt(CuArray, bccp_true)
        alpha_wage=CUDA.adapt(CuArray, alpha_wage)
        bwage=CUDA.adapt(CuArray, bwage)
        Lambda=CUDA.adapt(CuArray, Lambda)
        Pi_cpl=CUDA.adapt(CuArray,Pi_cpl)
        Pi_single_m=CUDA.adapt(CuArray,Pi_single_m)
        Pi_single_f=CUDA.adapt(CuArray,Pi_single_f)
        Pi_single=CUDA.adapt(CuArray,Pi_single)
        xval=CUDA.adapt(CuArray,xval)
        zval=CUDA.adapt(CuArray,zval)
        hval=CUDA.adapt(CuArray,hval)

        placeholders(b1, dt, dt_rf)
end

function f_convert_from_GPU()
    global dt, dt_rf, dt_rf_grid, Xstruct_grid, CCP_true, dt_struct_GMM, b1, PType, PType_true,
    xtran, xtran_married, bccp, bccp_true, Lambda, Pi_cpl, Pi_single, Pi_single_m, Pi_single_f,  xval, hval, 
    zval, alpha_wage, bwage, 
    dt_CPU, dt_rf_CPU, dt_rf_grid_CPU, Xstruct_grid_CPU, CCP_true_CPU, dt_struct_GMM_CPU, xval_CPU, zval_CPU, hval_CPU,
    b1_CPU, PType_CPU, PType_true_CPU, xtran_CPU, xtran_married_CPU, bccp_CPU, bccp_true_CPU, Lambda_CPU, Pi_cpl_CPU, Pi_single_m_CPU, Pi_single_f_CPU, alpha_wage_CPU, bwage_CPU, Pi_single_CPU

    global zeros_dim_x_st2_2, dummy_size0, dummy_size10, dummy_sizeNT, 
    zeros_st1_m_part1, zeros_st1_m_part2, zeros_st1_f_part1, zeros_st1_f_part2, zeros_st2_cpl_m, zeros_st2_cpl_f, zeros_st2_m, 
    zeros_st2_f, zeros_dim_x_st2_cpl_2, ones_st2_married, ones_st2_single

    dt= dt_CPU
    dt_rf=dt_rf_CPU
    dt_rf_grid=dt_rf_grid_CPU
    Xstruct_grid=Xstruct_grid_CPU
    CCP_true=CCP_true_CPU
    dt_struct_GMM=dt_struct_GMM_CPU
    xval=xval_CPU
    zval=zval_CPU
    hval=hval_CPU
    b1=b1_CPU
    PType=PType_CPU
    PType_true=PType_true_CPU
    xtran=xtran_CPU
    xtran_married=xtran_married_CPU
    bccp=bccp_CPU
    bccp_true=bccp_true_CPU
    alpha_wage=alpha_wage_CPU
    bwage=bwage_CPU
    Lambda=Lambda_CPU
    Pi_cpl=Pi_cpl_CPU
    Pi_single_m=Pi_single_m_CPU
    Pi_single_f=Pi_single_f_CPU
    Pi_single=Pi_single_CPU

    placeholders(b1, dt, dt_rf)
end