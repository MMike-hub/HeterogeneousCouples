
# This specification of the structural CCPs allows for gender-specific utility parameters
function CCP_structural(bccp_,Lambda, Xstruct_grid, fvt1_RX1, CCP_rf, i=0, bccp_dict_pseudoMLE=[])
    global dummy_size0

    p0_st1_m_marriedb4=dummy_size0
    p0_st1_f_marriedb4=dummy_size0
    p0_st1_m_singleb4=dummy_size0
    p0_st1_f_singleb4=dummy_size0
    p0_st2_m=dummy_size0
    p0_st2_f=dummy_size0
    p0_st2_cpl=dummy_size0
    p_st1_m_marriedb4=dummy_size0
    p_st1_f_marriedb4=dummy_size0
    p_st1_m_singleb4=dummy_size0
    p_st1_f_singleb4=dummy_size0
    p_st2_m=dummy_size0
    p_st2_f=dummy_size0
    p_st2_cpl=dummy_size0

    dim_Xstruct_st1_marriedb4_1=Xstruct_grid.dim_Xstruct_st1_marriedb4_1
    dim_Xstruct_st1_singleb4_1=Xstruct_grid.dim_Xstruct_st1_singleb4_1
    dim_Xstruct_st2_single_1= Xstruct_grid.dim_Xstruct_st2_single_1
    dim_Xstruct_st2_married_1=Xstruct_grid.dim_Xstruct_st2_married_1

    if i == 0
        bccp_st1_m_marriedb4 = get_slice(bccp_, 1, dim_Xstruct_st1_marriedb4_1)
        bccp_st1_m_singleb4 = bccp_st1_m_marriedb4
        bccp_st1_f_marriedb4 = get_slice(bccp_, dim_Xstruct_st1_marriedb4_1+1, dim_Xstruct_st1_marriedb4_1)
        bccp_st1_f_singleb4 = bccp_st1_f_marriedb4
        bccp_st2_m_single = get_slice(bccp_, 2*(dim_Xstruct_st1_marriedb4_1)+1, dim_Xstruct_st2_single_1)
        bccp_st2_f_single = get_slice(bccp_,  2*(dim_Xstruct_st1_marriedb4_1)+dim_Xstruct_st2_single_1+1, dim_Xstruct_st2_single_1)
        bccp_st2_m_married = get_slice(bccp_, 2*(dim_Xstruct_st1_marriedb4_1 +dim_Xstruct_st2_single_1)+1, dim_Xstruct_st2_married_1)
        bccp_st2_f_married = get_slice(bccp_, 2*(dim_Xstruct_st1_marriedb4_1 +dim_Xstruct_st2_single_1) + dim_Xstruct_st2_married_1+1, dim_Xstruct_st2_married_1)
    elseif i == 1
        bccp_st1_m_marriedb4 = get_slice(bccp_, 1, dim_Xstruct_st1_marriedb4_1)
        bccp_st1_m_singleb4 = bccp_st1_m_marriedb4
    elseif i == 2
        bccp_st1_f_marriedb4 = bccp_
        bccp_st1_f_singleb4 = bccp_st1_f_marriedb4
    elseif i == 3
        bccp_st2_m_single = bccp_
    elseif i == 4
        bccp_st2_f_single = bccp_
    elseif i == 5
        bccp_st2_m_married = get_slice(bccp_, 1, dim_Xstruct_st2_married_1)
        bccp_st2_f_married = get_slice(bccp_, dim_Xstruct_st2_married_1+1, dim_Xstruct_st2_married_1)
    end
  
    # First-stage CCPs depend on second stage married utilities. So when I am estimating the first-stage parameters, I pass the relevant second-stage parameters as fixed parameters.
    # In order to ensure that at every EM loop, first-stage CCPs are produced using the most updated version of the second-stage CCPs, I make sure that MLE over the second-stage parameters runs before the MLE over first-stage parameters
    if (i==1) | (i==5) | (i==0)
        if (i==1)
            bccp_st2_m_married=get_slice(bccp_dict_pseudoMLE[5], 1, dim_Xstruct_st2_married_1)
        end
        du_m_married=sum(Xstruct_grid.Xstruct_st2_married .* bccp_st2_m_married,dims=1)
        du_m_married=reshape(du_m_married,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))
    end
    if (i==2) | (i==5) | (i==0)
        if (i==2)
            bccp_st2_f_married=get_slice(bccp_dict_pseudoMLE[5], dim_Xstruct_st2_married_1+1, dim_Xstruct_st2_married_1)
        end
        du_f_married=sum(Xstruct_grid.Xstruct_st2_married .* bccp_st2_f_married,dims=1)
        du_f_married=reshape(du_f_married,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))
    end
    if (i==1) | (i==0)
        # First stage structural CCPs.
        dkappa_m_marriedb4= sum(Xstruct_grid.Xstruct_st1_marriedb4 .* bccp_st1_m_marriedb4,dims=1) # This is the DIFFERENCE in kappas with respect to remaining single, of course, I cannot identify the levels of kappa
        dkappa_m_marriedb4=reshape(dkappa_m_marriedb4,(phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin, hbin, NTypes, phibin_sp_orig, hbin_sp_orig, NTypes_sp_orig, T2))

        # Generate that part of the continuation value that depends on utilities
        du_m_married_fv = sum(CCP_rf.p_st2_cpl .* du_m_married, dims=(7,8))
        # Drop singleton dimensions
        du_m_married_fv= dropdims(du_m_married_fv, dims=(7,8))
        du_m_married_fv_marriedb4_view = @view du_m_married_fv[Xstruct_grid.idx_st2_into_st1];
        du_m_married_fv_marriedb4_view= reshape(du_m_married_fv_marriedb4_view, (phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin_m, hbin_m, NTypes_m, 1,1,1, T2))
        dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[Xstruct_grid.idx_st2_into_st1]
        dEU_st1_m_view = reshape(dEU_st1_m_view, (phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin_m, hbin_m, NTypes_m, 1,1,1, T2))
        # du_m_married_fv_marriedb4_view= @view du_m_married_fv[phistate_list[:,1],:,:,phistate_list[:,1],:,:,:,:,:,:] # Having two phistate_list[:,1] indexes in different fields doesn't work. The view is created, but then I cannot perform any operations with it (when du_m_married_fv is CuArray)
        # dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[phistate_list[:,1],:,:,phistate_list[:,1],:,:,:] # Same here
        # I have to add some placeholder dimensions to dEU_st1_m, those corresponding to the origin spouse
        num_p_st1_m_marriedb4= exp.(dkappa_m_marriedb4  .+ du_m_married_fv_marriedb4_view .+ dEU_st1_m_view)
        den_p_st1_m_marriedb4=sum(num_p_st1_m_marriedb4,dims=(1,2,3)) .+ 1 
        p_st1_m_marriedb4=num_p_st1_m_marriedb4./den_p_st1_m_marriedb4
        p0_st1_m_marriedb4 = dropdims(1 ./den_p_st1_m_marriedb4,dims=(1,2,3))

        dkappa_m_singleb4= sum(Xstruct_grid.Xstruct_st1_singleb4 .* bccp_st1_m_singleb4,dims=1)
        dkappa_m_singleb4=reshape(dkappa_m_singleb4,(phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin, hbin, NTypes, T2))

        du_m_married_fv_singleb4_view= @view du_m_married_fv[Xstruct_grid.idx_st2_into_st1]
        dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[Xstruct_grid.idx_st2_into_st1]
        num_p_st1_m_singleb4=exp.(dkappa_m_singleb4 .+ du_m_married_fv_singleb4_view  .+  dEU_st1_m_view)
        den_p_st1_m_singleb4=sum(num_p_st1_m_singleb4,dims=(1,2,3)) .+ 1
        p_st1_m_singleb4=num_p_st1_m_singleb4./den_p_st1_m_singleb4
        p0_st1_m_singleb4= dropdims(1 ./den_p_st1_m_singleb4,dims=(1,2,3))
    end

    if (i==2) | (i==0)
        dkappa_f_marriedb4= sum(Xstruct_grid.Xstruct_st1_marriedb4 .* bccp_st1_f_marriedb4,dims=1) # This is the DIFFERENCE in kappas with respect to remaining single, of course, I cannot identify the levels of kappa
        dkappa_f_marriedb4=reshape(dkappa_f_marriedb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))
        du_f_married_fv = sum(CCP_rf.p_st2_cpl .* du_f_married, dims=(7,8))
        du_f_married_fv= dropdims(du_f_married_fv, dims=(7,8))

        du_f_married_fv_marriedb4_view = @view du_f_married_fv[Xstruct_grid.idx_st2_into_st1];
        du_f_married_fv_marriedb4_view= reshape(du_f_married_fv_marriedb4_view, (phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,1,1,1,T2))
        dEU_st1_f_view= @view fvt1_RX1.dEU_st1_m[Xstruct_grid.idx_st2_into_st1]
        dEU_st1_f_view = reshape(dEU_st1_f_view, (phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,1,1,1,T2))
        num_p_st1_f_marriedb4= exp.(dkappa_f_marriedb4 .+ du_f_married_fv_marriedb4_view  .+ dEU_st1_f_view)
        den_p_st1_f_marriedb4=sum(num_p_st1_f_marriedb4,dims=(1,2,3)) .+ 1
        p_st1_f_marriedb4=num_p_st1_f_marriedb4./den_p_st1_f_marriedb4
        p0_st1_f_marriedb4= dropdims(1 ./ den_p_st1_f_marriedb4,dims=(1,2,3))

        dkappa_f_singleb4= sum(Xstruct_grid.Xstruct_st1_singleb4 .* bccp_st1_f_singleb4,dims=1)
        dkappa_f_singleb4=reshape(dkappa_f_singleb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))
        du_f_married_fv_singleb4_view= @view du_f_married_fv[Xstruct_grid.idx_st2_into_st1]
        dEU_st1_f_view= @view fvt1_RX1.dEU_st1_f[Xstruct_grid.idx_st2_into_st1]
        num_p_st1_f_singleb4=exp.(dkappa_f_singleb4 .+ du_f_married_fv_singleb4_view .+ dEU_st1_f_view)
        den_p_st1_f_singleb4=sum(num_p_st1_f_singleb4,dims=(1,2,3)) .+ 1
        p_st1_f_singleb4=num_p_st1_f_singleb4./den_p_st1_f_singleb4
        p0_st1_f_singleb4=  dropdims(1 ./den_p_st1_f_singleb4,dims=(1,2,3))
    end

    if (i==3) | (i==0)
        # Second stage structural CCPs
        du_m_single=sum(Xstruct_grid.Xstruct_st2_single .* bccp_st2_m_single,dims=1)
        du_m_single=reshape(du_m_single,(xbin_m,hbin_m,NTypes_m,J_m+1,T2))
        num_p_st2_m_single=exp.(du_m_single   .+ fvt1_RX1.dEU_st2_single_m)
        den_p_st2_m_single=sum(num_p_st2_m_single,dims=4) 
        p_st2_m=num_p_st2_m_single./den_p_st2_m_single
        p0_st2_m=p_st2_m[:,:,:,1,:]
    end
    if (i==4) | (i==0)
        du_f_single=sum(Xstruct_grid.Xstruct_st2_single .* bccp_st2_f_single,dims=1)
        du_f_single=reshape(du_f_single,(xbin_m,hbin_m,NTypes_m,J_m+1,T2))
        num_p_st2_f_single=exp.(du_f_single   .+ fvt1_RX1.dEU_st2_single_m)
        den_p_st2_f_single=sum(num_p_st2_f_single,dims=4) 
        p_st2_f=num_p_st2_f_single./den_p_st2_f_single
        p0_st2_f=p_st2_f[:,:,:,1,:]
    end
    if (i==5) | (i==0)
        dv_m_married=du_m_married .+ fvt1_RX1.dEU_st2_married_m
        dv_f_married=du_f_married .+ fvt1_RX1.dEU_st2_married_f

        num_p_st2_cpl= exp.(Lambda .* dv_m_married .+ (1 .- Lambda) .* dv_f_married)
        den_p_st2_cpl= sum(num_p_st2_cpl,dims=(7,8))
        p_st2_cpl=num_p_st2_cpl ./ den_p_st2_cpl
        p0_st2_cpl=p_st2_cpl[:,:,:,:,:,:,1,1,:]
    end


    CCP_struct= CCP_Structural(p0_st1_m_marriedb4, 
    p0_st1_f_marriedb4, 
    p0_st1_m_singleb4, 
    p0_st1_f_singleb4, 
    p0_st2_m, 
    p0_st2_f, 
    p0_st2_cpl, 
    p_st1_m_marriedb4, 
    p_st1_f_marriedb4, 
    p_st1_m_singleb4, 
    p_st1_f_singleb4, 
    p_st2_m, 
    p_st2_f, 
    p_st2_cpl)


    return CCP_struct
end


struct CCP_Structural
    p0_st1_m_marriedb4:: Union{Array, CuArray} 
    p0_st1_f_marriedb4:: Union{Array, CuArray} 
    p0_st1_m_singleb4:: Union{Array, CuArray} 
    p0_st1_f_singleb4:: Union{Array, CuArray} 
    p0_st2_m:: Union{Array, CuArray} 
    p0_st2_f:: Union{Array, CuArray} 
    p0_st2_cpl:: Union{Array, CuArray} 
    p_st1_m_marriedb4:: Union{Array, CuArray} 
    p_st1_f_marriedb4:: Union{Array, CuArray} 
    p_st1_m_singleb4:: Union{Array, CuArray} 
    p_st1_f_singleb4:: Union{Array, CuArray} 
    p_st2_m:: Union{Array, CuArray} 
    p_st2_f:: Union{Array, CuArray} 
    p_st2_cpl:: Union{Array, CuArray} 
end

# Use views to avoid copying slices
function get_slice(data, start, len)
    return @view data[start:(start+len-1)]
end