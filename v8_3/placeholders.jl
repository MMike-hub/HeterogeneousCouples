function placeholders(bccp, dt)
    global zeros_dim_x_st2_2, dummy_size0, dummy_size10, dummy_sizeNT, 
    zeros_st1_m_part1, zeros_st1_m_part2, zeros_st1_f_part1, zeros_st1_f_part2, zeros_st2_cpl_m, zeros_st2_cpl_f, zeros_st2_m, 
    zeros_st2_f, zeros_dim_x_st2_cpl_2
    # Allocate memory for placeholders. Can't do it inside the objective function or Zygote will throw an error
    if bccp isa CuArray
        
        # PLACEHOLDERS USED IN LikeCCP
        # Dimensionless array placeholder
        dummy_size0=CuArray{precision}(zeros(0))
        # Array with same size as data placeholder
        dummy_sizeNT=CUDA.fill(convert(precision, 1), size(dt.idx_p_st1_m_singleb4))
        # Columnless matrix placeholder
        dummy_size10=CUDA.fill(convert(precision, 0), (1,0))
        
        #PLACEHOLDERS USED IN fvdata_RX1
        zeros_st1_m_part1=CUDA.fill(convert(precision, 0), (xbin_m_t,xbin_sp_t,hbin_m,hbin_sp,NTypes_m,NTypes_sp_orig,J_m+1,J_sp+1,1))
        zeros_st1_m_part2=CUDA.fill(convert(precision, 0), (xbin_m,hbin_m,NTypes_m,1))
        zeros_st1_f_part1 = CUDA.fill(convert(precision, 0), (xbin_sp_t,xbin_f_t,hbin_sp,hbin_f,NTypes_sp_orig,NTypes_f,J_sp+1,J_f+1,1))
        zeros_st1_f_part2 = CUDA.fill(convert(precision, 0), (xbin_f,hbin_f,NTypes_f,1))
        zeros_st2_cpl_m= CUDA.fill(convert(precision, 0), (xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_sp_orig,J_m+1,J_f+1,1))
        zeros_st2_cpl_f= CUDA.fill(convert(precision, 0), (xbin_m,xbin_f,hbin_m,hbin_f,NTypes_sp_orig,NTypes_f,J_m+1,J_f+1,1))
        zeros_st2_m= CUDA.fill(convert(precision, 0), (xbin_m,hbin_m,NTypes_m,J_m+1,1))
        zeros_st2_f= CUDA.fill(convert(precision, 0), (xbin_f,hbin_f,NTypes_f,J_f+1,1))
    else
        
        # PLACEHOLDERS USED IN LikeCCP
        # Dimensionless array placeholder
        dummy_size0=zeros(precision, 0)
        # Array with same size as data placeholder
        dummy_sizeNT=ones(precision, size(dt.idx_p_st1_m_singleb4))
        dummy_size10=zeros(precision , (1,0))
        
        zeros_st1_m_part1=zeros(precision, (xbin_m_t,xbin_sp_t,hbin_m,hbin_sp,NTypes_m,NTypes_sp_orig,J_m+1,J_sp+1,1))
        zeros_st1_m_part2=zeros(precision, (xbin_m,hbin_m,NTypes_m,1))
        zeros_st1_f_part1 = zeros(precision, (xbin_sp_t,xbin_f_t,hbin_sp,hbin_f,NTypes_sp_orig,NTypes_f,J_sp+1,J_f+1,1))
        zeros_st1_f_part2 = zeros(precision, (xbin_f,hbin_f,NTypes_f,1))
        zeros_st2_cpl_m= zeros(precision, (xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_sp_orig,J_m+1,J_f+1,1))
        zeros_st2_cpl_f= zeros(precision, (xbin_m,xbin_f,hbin_m,hbin_f,NTypes_sp_orig,NTypes_f,J_m+1,J_f+1,1))
        zeros_st2_m= zeros(precision, (xbin_m,hbin_m,NTypes_m,J_m+1,1))
        zeros_st2_f= zeros(precision, (xbin_f,hbin_f,NTypes_f,J_f+1,1))
    end
end