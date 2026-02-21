function gmm_struct(bccp, dt, Lambda, Xstruct_grid, fvt1_RX1, PType, W, i=1)
    CCP_struct=CCP_structural(bccp, Lambda, Xstruct_grid, fvt1_RX1)
    res_2=resCCP(dt,CCP_struct)

    Z_st1_m_singleb4_data= reshape([dt.x dt.z dt.s zeros(size(dt.x)) zeros(size(dt.x)) zeros(size(dt.x))],(:,1,6))
    prod_res_st1_m_singleb4_data= reshape(Z_st1_m_singleb4_data .* res.res_st1_m_singleb4_data, (length(dt.x),:))
    prod_res_st1_m_singleb4_data=coalesce.(prod_res_st1_m_singleb4_data,0) .* (dt.sex.===1) 

    Z_st1_f_singleb4_data= reshape([dt.x dt.z dt.s zeros(size(dt.x)) zeros(size(dt.x)) zeros(size(dt.x))],(:,1,6))
    prod_res_st1_f_singleb4_data= reshape(Z_st1_f_singleb4_data .* res.res_st1_f_singleb4_data, (length(dt.x),:))
    prod_res_st1_f_singleb4_data=coalesce.(prod_res_st1_f_singleb4_data,0) .* (dt.sex.===2)
    
    Z_st1_m_marriedb4_data= reshape([dt.x dt.z dt.s dt.x_sp_orig dt.z_sp_orig dt.s_sp_orig],(:,1,6))
    prod_res_st1_m_marriedb4_data= reshape(Z_st1_m_marriedb4_data .* res.res_st1_m_marriedb4_data, (length(dt.x),:))
    prod_res_st1_m_marriedb4_data=coalesce.(prod_res_st1_m_marriedb4_data,0) .* (dt.sex.===1)
    Z_st1_f_marriedb4_data= reshape([dt.x dt.z dt.s dt.x_sp_orig dt.z_sp_orig dt.s_sp_orig],(:,1,6))
    prod_res_st1_f_marriedb4_data= reshape(Z_st1_f_marriedb4_data .* res.res_st1_f_marriedb4_data, (length(dt.x),:))
    prod_res_st1_f_marriedb4_data=coalesce.(prod_res_st1_f_marriedb4_data,0) .* (dt.sex.===2)

    prod_res_st1_f_data=prod_res_st1_f_marriedb4_data .* coalesce.(dt.married_tminus1,0) .+ prod_res_st1_f_singleb4_data .*(1 .- coalesce.(dt.married_tminus1,0))
    prod_res_st1_m_data=prod_res_st1_m_marriedb4_data .* coalesce.(dt.married_tminus1,0) .+ prod_res_st1_m_singleb4_data .*(1 .- coalesce.(dt.married_tminus1,0))

    prod_res_st1_data=prod_res_st1_f_data .* (dt.sex.===2) .+ prod_res_st1_m_data .* (dt.sex.===1)
    
    Z_st2_f_data= reshape([dt.x dt.z dt.s],(:,1,3))
    prod_res_st2_f_data=reshape(Z_st2_f_data .* res.res_st2_f_data, (length(dt.x),:))
    prod_res_st2_f_data=coalesce.(prod_res_st2_f_data,0).*(dt.sex.===2)

    Z_st2_m_data= reshape([dt.x dt.z dt.s],(:,1,3))
    prod_res_st2_m_data=reshape(Z_st2_m_data .* res.res_st2_m_data, (length(dt.x),:))
    prod_res_st2_m_data=coalesce.(prod_res_st2_m_data,0).*(dt.sex.===1)

    prod_res_st2_single_data=(prod_res_st2_f_data .* (dt.sex.===2) .+ prod_res_st2_m_data .* (dt.sex.===1)).*(dt.married.===0)

    Z_st2_cpl_data=reshape([dt.x_m_cpl dt.z_m_cpl dt.s_m_cpl dt.x_f_cpl dt.z_f_cpl dt.s_f_cpl],(:,1,6))
    prod_res_st2_cpl_data=reshape(Z_st2_cpl_data .* res.res_st2_cpl_data, (length(dt.x_m_cpl),:))
    prod_res_st2_cpl_data=coalesce.(prod_res_st2_cpl_data,0).*(dt.married.===1)

    # I can assume that the moment condition is the orthogonality between each instrument and 
    # the residual corresponding to each choice, minus one, i.e. E(Z*e_j)=0 for all choices j=1,2 (where j=0 is engine replacement)
    # Where e_j := d_j - p_j(x)
    # This means that there are as many moment conditions as there are choices (minus one) times the number 
    # of state variables whose parameter we need to estimate. That's J*dim_X.
    # The problem with this is that the covariance matrix becomes singular.
    prod_ =[prod_res_st1_data prod_res_st2_single_data prod_res_st2_cpl_data]

    # Instead, I can show that the FOC of the MLE estimator for a multinomial logit where the parameters are 
    # identical across choices and it's the regressors that change across choices, corresponds to the moment condition
    # ∑_j=1...J E[X_dis{jk}(1(d_{ntj}==1) - p_j(X))]=0 ,  where X_dis{jk} is the value of the state variable that multiplies
    # the parameter β_k in the utility corresponding to choice j. SEE MISC NOTES 3 search for "Flexible logit MLE FOC"
    # So, all I have to do is to sum prod_ *across choices*
    # prod_=sum(reshape(prod_,(dim, dim_X, J)),dims=3)[:,:,1]
    prod_weighted=PType .* prod_

    # Sum across types
    dim_prod1=size(dt.x,1)
    dim_prod2=size(prod_weighted,2)
    dim_s=Int(dim_prod1/NTypes)
    prod_weighted=reshape(prod_weighted,(dim_s,NTypes,dim_prod2))
    prod_weighted=permutedims(prod_weighted,(1,3,2))
    prod_weighted=sum(prod_weighted,dims=3)[:,:,1]

    denom_prod_st1_data=length(dt.x)*ones(size(prod_res_st1_data,2))
    denom_prod_st2_single_data=sum(coalesce.(dt.married,0).==0)*ones(size(prod_res_st2_single_data,2))
    denom_prod_st2_cpl_data=sum(coalesce.(dt.married,0).==1)*ones(size(prod_res_st2_cpl_data,2))

    denom_prod=vcat(denom_prod_st1_data , denom_prod_st2_single_data, denom_prod_st2_cpl_data)

    if i in [1,  3]
        moment=sum(prod_weighted,dims=1)[:]./ denom_prod
        if i==1
            obj=sum(moment' * W * moment)
        else
            obj=[]
        end
    else
        moment=[]
        obj=[]
    end

    if i ==2
        Σ=prod_weighted' * prod_weighted ./ denom_prod
    else
        Σ=[]
    end    

    return obj, Σ, moment, prod_weighted
end