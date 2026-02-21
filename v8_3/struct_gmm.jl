function f_struct_gmm(bccp_, dt_struct_GMM, Xstruct_grid, fvt1_RX1, CCP_, PType, W, i, j=1)

    #For testing
    # b1_=ones(dt_struct_GMM.dim_z_st1_2*dt_struct_GMM.dim_y_st1+dt_struct_GMM.dim_z_st1_2*dt_struct_GMM.dim_y_st1+dt_struct_GMM.dim_z_st2_2*dt_struct_GMM.dim_y_st2+dt_struct_GMM.dim_z_st2_2*dt_struct_GMM.dim_y_st2+dt_struct_GMM.dim_z_st2_cpl_2*dt_struct_GMM.dim_y_st2_cpl)
    
    CCP_struct=CCP_structural(bccp_ ,Lambda, Xstruct_grid, fvt1_RX1, CCP_, i, bccp_dict_pseudoMLE)
    res=resCCP(dt, CCP_struct, i)

    dim=dt_struct_GMM.dim_z_1
    if [0,1,3]
        # Z_st1_m_singleb4_data= reshape(dt_struct_GMM.Z_st1_m_singleb4_data,(:,1,dt_struct_GMM.dim_z_st1_2))
        # prod_res_st1_m_singleb4_data= reshape(Z_st1_m_singleb4_data .* res.res_st1_m_singleb4_data, (dim,:))
        # prod_res_st1_m_singleb4_data=coalesce.(prod_res_st1_m_singleb4_data,0) .* (dt_struct_GMM.sex.===1) .* (1 .- dt.firstsampled)
        
        # Z_st1_m_marriedb4_data= reshape(dt_struct_GMM.Z_st1_m_marriedb4_data,(:,1,dt_struct_GMM.dim_z_st1_2))
        # prod_res_st1_m_marriedb4_data= reshape(Z_st1_m_marriedb4_data .* res.res_st1_m_marriedb4_data, (dim,:))
        # prod_res_st1_m_marriedb4_data=coalesce.(prod_res_st1_m_marriedb4_data,0) .* (dt_struct_GMM.sex.===1) .* (1 .- dt.firstsampled)

        # prod_res_st1_m_data=similar(prod_res_st1_m_marriedb4_data)
        # prod_res_st1_m_data.=prod_res_st1_m_marriedb4_data .* dt.married_tminus1_ .+ prod_res_st1_m_singleb4_data .*(1 .- dt.married_tminus1_)
        
        Z_st1_m_data=reshape(dt_struct_GMM.Z_st1_m_singleb4_data,(:,1,dt_struct_GMM.dim_z_st1_2))
        prod_res_st1_m_data=(Z_st1_m_data .* res.res_st1_m_marriedb4_data ) .* dt.married_tminus1_ .+ (Z_st1_m_data .* res.res_st1_m_singleb4_data) .*(1 .- dt.married_tminus1_)
        prod_res_st1_m_data= prod_res_st1_m_data .*= ((dt.sex_.==1) .* (1 .- dt.firstsampled))
        prod_res_st1_m_data=reshape(prod_res_st1_m_data, (dim,:))

        denom_prod_st1_m_data=sum(dt.sex_.==1)*ones(size(prod_res_st1_m_data,2))    
    end
    if [0,2,3]
        Z_st1_f_singleb4_data= reshape(dt_struct_GMM.Z_st1_f_singleb4_data,(:,1,dt_struct_GMM.dim_z_st1_2))
        prod_res_st1_f_singleb4_data= reshape(Z_st1_f_singleb4_data .* res.res_st1_f_singleb4_data, (dim,:))
        prod_res_st1_f_singleb4_data=coalesce.(prod_res_st1_f_singleb4_data,0) .* (dt_struct_GMM.sex.===2) .* (1 .- dt.firstsampled)
    
        Z_st1_f_marriedb4_data= reshape(dt_struct_GMM.Z_st1_f_marriedb4_data,(:,1,dt_struct_GMM.dim_z_st1_2))
        prod_res_st1_f_marriedb4_data= reshape(Z_st1_f_marriedb4_data .* res.res_st1_f_marriedb4_data, (dim,:))
        prod_res_st1_f_marriedb4_data=coalesce.(prod_res_st1_f_marriedb4_data,0) .* (dt_struct_GMM.sex.===2) .* (1 .- dt.firstsampled)
        
        prod_res_st1_f_data=similar(prod_res_st1_f_marriedb4_data)
        prod_res_st1_f_data=prod_res_st1_f_marriedb4_data .* dt.married_tminus1_ .+ prod_res_st1_f_singleb4_data .*(1 .- dt.married_tminus1_)
       
        denom_prod_st1_f_data=sum(dt.sex_.==2)*ones(size(prod_res_st1_f_data,2))
        denom_prod=denom_prod_st1_f_data
        
        prod_weighted=similar(prod_res_st1_f_data)
        prod_weighted =PType .*prod_res_st1_f_data.*dt.insample
    end
    if [0,4]
        Z_st2_m_data=reshape(dt_struct_GMM.Z_st2_m_data,(:,1,dt_struct_GMM.dim_z_st2_2))
        prod_res_st2_m_data=reshape(Z_st2_m_data .* res.res_st2_m_data, (dim,:))
        prod_res_st2_m_data=coalesce.(prod_res_st2_m_data,0).*(dt_struct_GMM.sex.===1)
        denom_prod_st2_m_data=sum((dt.married_.==0).*(dt.sex_.==1))*ones(size(prod_res_st2_m_data,2))
        denom_prod=denom_prod_st2_m_data
        prod_weighted =PType .*prod_res_st2_m_data.*dt.insample

        Z_st2_f_data= reshape(dt_struct_GMM.Z_st2_f_data,(:,1,dt_struct_GMM.dim_z_st2_2))
        prod_res_st2_f_data=reshape(Z_st2_f_data .* res.res_st2_f_data, (dim,:))
        prod_res_st2_f_data=coalesce.(prod_res_st2_f_data,0).*(dt_struct_GMM.sex.===2)
        denom_prod_st2_f_data=sum((dt.married_.==0).*(dt.sex_.==2))*ones(size(prod_res_st2_f_data,2))
        denom_prod=denom_prod_st2_f_data
        prod_weighted =PType .*prod_res_st2_f_data.*dt.insample

        Z_st2_cpl_data=reshape(dt_struct_GMM.Z_st2_cpl_data,(:,1,dt_struct_GMM.dim_z_st2_cpl_2))
        prod_res_st2_cpl_data=reshape(Z_st2_cpl_data .* res.res_st2_cpl_data, (dim,:))
        prod_res_st2_cpl_data=coalesce.(prod_res_st2_cpl_data,0).*(dt_struct_GMM.married.===1)
        denom_prod_st2_cpl_data=sum(dt.married_.==1)*ones(size(prod_res_st2_cpl_data,2))
        denom_prod=denom_prod_st2_cpl_data
        prod_weighted =PType .*prod_res_st2_cpl_data.*dt.insample 
    end
    

    # Sum across types
    dim_prod2=size(prod_weighted,2)
    dim_s=Int(dim/NTypes)
    prod_weighted=reshape(prod_weighted,(dim_s,NTypes,dim_prod2))
    prod_weighted=permutedims(prod_weighted,(1,3,2))
    prod_weighted=sum(prod_weighted,dims=3)[:,:,1]


    if j in [1, 3]
        moment_dis_rf=sum(prod_weighted,dims=1)[:] ./ denom_prod
        if j ==1
            # sum() only used to convert 1 by 1 matrix into scalar
            obj=sum(moment_dis_rf' * W * moment_dis_rf)
        else 
            obj=[]
        end
    else
        moment_dis_rf=[]
        obj=[]
    end
    if j==2
        Σ=prod_weighted' * prod_weighted ./ denom_prod
    else
        Σ=[]
    end
    
    return obj, Σ, moment_dis_rf, prod_weighted
end