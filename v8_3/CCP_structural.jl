# Unlike the case of the reduced form regression, 
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
    dim_Xstruct_st2_m_married_1=Xstruct_grid.dim_Xstruct_st2_m_married_1
    dim_Xstruct_st2_f_married_1=Xstruct_grid.dim_Xstruct_st2_f_married_1

    # During estimation, it will never be useful to calculate 2nd stage CCPs for men, women, and couples separately, since they all depend on the same parameters. 
    # i==0 --> simultaneously calculate CCPs for of all choices
    # i==1 --> only calculate CCPs for male's 1st-stage choices
    # i==2 --> only calculate CCPs for female's 1st-stage choices
    # i==3 --> only calculate CCPs for male and female's 1st-stage choices
    # i==4 --> only calculate CCPs for all 2nd-stage choices

    if (i==0)
        bccp_st1_m_marriedb4, bccp_st1_m_singleb4, bccp_st1_f_marriedb4, bccp_st1_f_singleb4, bccp_st2_m_single, bccp_st2_f_single, bccp_st2_m_married, bccp_st2_f_married= bccp_split(bccp_, Xstruct_grid)
    elseif (i==1)
        bccp_st1_m_marriedb4 = bccp_
        bccp_st1_m_singleb4 = bccp_st1_m_marriedb4
    elseif (i==2)
        bccp_st1_f_marriedb4 = bccp_
        bccp_st1_f_singleb4 = bccp_st1_f_marriedb4
    elseif (i==3)
        bccp_st1_m_marriedb4=bccp_
        bccp_st1_m_singleb4 = bccp_st1_m_marriedb4
        bccp_st1_f_marriedb4=bccp_
        bccp_st1_f_singleb4 = bccp_st1_f_marriedb4
    elseif (i==4)
        bccp_st2_m_single = bccp_
        bccp_st2_f_single = bccp_st2_m_single
        bccp_st2_m_married = bccp_st2_m_single
        bccp_st2_f_married = bccp_st2_m_single
    end
  
    # First-stage CCPs depend on second stage married utilities. So when I am estimating the first-stage parameters, I pass the relevant second-stage parameters as fixed parameters.
    # In order to ensure that at every EM loop, first-stage CCPs are produced using the most updated version of the second-stage CCPs, I make sure that MLE over the second-stage parameters runs before the MLE over first-stage parameters, that's why the order of the loop in wlogit_struct_wrap matters.
    if i in [0, 1, 3, 4] 
        if i in [1, 3]
            bccp_st2_m_married=bccp_dict_pseudoMLE[4]
        end
        du_m_married=sum(Xstruct_grid.Xstruct_st2_m_married .* bccp_st2_m_married,dims=1)
        du_m_married=reshape(du_m_married,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))
    end
    if i in [0, 2, 3, 4] 
        if i in [2, 3]
            bccp_st2_f_married=bccp_dict_pseudoMLE[4]
        end
        du_f_married=sum(Xstruct_grid.Xstruct_st2_f_married .* bccp_st2_f_married,dims=1)
        du_f_married=reshape(du_f_married,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))
    end
    if i in [0, 1, 3] 
        # First stage male structural CCPs.
        dkappa_m_marriedb4= sum(Xstruct_grid.Xstruct_st1_marriedb4 .* bccp_st1_m_marriedb4,dims=1) # This is the DIFFERENCE in kappas with respect to remaining single, of course, I cannot identify the levels of kappa
        dkappa_m_marriedb4=reshape(dkappa_m_marriedb4,(phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin, hbin, NTypes, phibin_sp_orig, hbin_sp_orig, NTypes_sp_orig, T2))

        # Generate that part of the continuation value that depends on utilities
        du_m_married_fv = sum(CCP_rf.p_st2_cpl .* du_m_married, dims=(7,8))
        # Drop singleton dimensions
        du_m_married_fv= dropdims(du_m_married_fv, dims=(7,8))
        du_m_married_fv_marriedb4_view = @view du_m_married_fv[Xstruct_grid.idx_st2_into_st1_m]
        du_m_married_fv_marriedb4_view= reshape(du_m_married_fv_marriedb4_view, (phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin_m, hbin_m, NTypes_m, 1,1,1, T2))
        dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[Xstruct_grid.idx_st2_into_st1_m] 
        dEU_st1_m_view = reshape(dEU_st1_m_view, (phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin_m, hbin_m, NTypes_m, 1,1,1, T2))
        # du_m_married_fv_marriedb4_view= @view du_m_married_fv[phistate_list[:,1],:,:,phistate_list[:,1],:,:,:,:,:,:] # Having two phistate_list[:,1] indexes in different fields doesn't work. The view is created, but then I cannot perform any operations with it (when du_m_married_fv is CuArray)
        # dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[phistate_list[:,1],:,:,phistate_list[:,1],:,:,:] # Same here
        # I have to add some placeholder dimensions to dEU_st1_m, those corresponding to the origin spouse
        max_=maximum(dkappa_m_marriedb4  .+ du_m_married_fv_marriedb4_view .+ dEU_st1_m_view)
        num_p_st1_m_marriedb4= exp.(dkappa_m_marriedb4  .+ du_m_married_fv_marriedb4_view .+ dEU_st1_m_view .- max_)
        den_p_st1_m_marriedb4= sum(num_p_st1_m_marriedb4,dims=(1,2,3)) .+ 1 
        p_st1_m_marriedb4= num_p_st1_m_marriedb4./den_p_st1_m_marriedb4 
        # maximum(abs.(p_st1_m_marriedb4.-CCP_true.p_st1_m_marriedb4)./CCP_true.p_st1_m_marriedb4)
        p0_st1_m_marriedb4 = dropdims(1 ./den_p_st1_m_marriedb4,dims=(1,2,3))

        dkappa_m_singleb4= sum(Xstruct_grid.Xstruct_st1_singleb4 .* bccp_st1_m_singleb4,dims=1)
        dkappa_m_singleb4=reshape(dkappa_m_singleb4,(phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin, hbin, NTypes, T2))

        du_m_married_fv_singleb4_view= @view du_m_married_fv[Xstruct_grid.idx_st2_into_st1_m]
        dEU_st1_m_view= @view fvt1_RX1.dEU_st1_m[Xstruct_grid.idx_st2_into_st1_m] #I know it's weird to use 
        max_=maximum(dkappa_m_singleb4 .+ du_m_married_fv_singleb4_view  .+  dEU_st1_m_view)
        num_p_st1_m_singleb4=exp.(dkappa_m_singleb4 .+ du_m_married_fv_singleb4_view  .+  dEU_st1_m_view .- max_)
        den_p_st1_m_singleb4=sum(num_p_st1_m_singleb4,dims=(1,2,3)) .+ 1
        p_st1_m_singleb4=num_p_st1_m_singleb4./den_p_st1_m_singleb4
        # maximum(abs.(p_st1_m_singleb4.-CCP_true.p_st1_m_singleb4)./CCP_true.p_st1_m_singleb4)
        p0_st1_m_singleb4= dropdims(1 ./den_p_st1_m_singleb4,dims=(1,2,3))
    end

    if i in [0, 2, 3] 
        # First stage female structural CCPs.
        dkappa_f_marriedb4= sum(Xstruct_grid.Xstruct_st1_marriedb4 .* bccp_st1_f_marriedb4,dims=1) # This is the DIFFERENCE in kappas with respect to remaining single, of course, I cannot identify the levels of kappa
        dkappa_f_marriedb4=reshape(dkappa_f_marriedb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))
        du_f_married_fv = sum(CCP_rf.p_st2_cpl .* du_f_married, dims=(7,8))
        du_f_married_fv= dropdims(du_f_married_fv, dims=(7,8))

        du_f_married_fv_marriedb4_view = @view du_f_married_fv[Xstruct_grid.idx_st2_into_st1_f];
        du_f_married_fv_marriedb4_view= reshape(du_f_married_fv_marriedb4_view, (phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,1,1,1,T2))
        dEU_st1_f_view= @view fvt1_RX1.dEU_st1_f[Xstruct_grid.idx_st2_into_st1_f]
        dEU_st1_f_view = reshape(dEU_st1_f_view, (phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,1,1,1,T2))
        max_=maximum(dkappa_f_marriedb4 .+ du_f_married_fv_marriedb4_view  .+ dEU_st1_f_view)
        num_p_st1_f_marriedb4= exp.(dkappa_f_marriedb4 .+ du_f_married_fv_marriedb4_view  .+ dEU_st1_f_view .- max_)
        den_p_st1_f_marriedb4=sum(num_p_st1_f_marriedb4,dims=(1,2,3)) .+ 1
        p_st1_f_marriedb4=num_p_st1_f_marriedb4./den_p_st1_f_marriedb4
        # maximum(abs.(p_st1_f_marriedb4.-CCP_true.p_st1_f_marriedb4)./CCP_true.p_st1_f_marriedb4)
        p0_st1_f_marriedb4= dropdims(1 ./ den_p_st1_f_marriedb4,dims=(1,2,3))

        dkappa_f_singleb4= sum(Xstruct_grid.Xstruct_st1_singleb4 .* bccp_st1_f_singleb4,dims=1)
        dkappa_f_singleb4=reshape(dkappa_f_singleb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))
        du_f_married_fv_singleb4_view= @view du_f_married_fv[Xstruct_grid.idx_st2_into_st1_f]
        dEU_st1_f_view= @view fvt1_RX1.dEU_st1_f[Xstruct_grid.idx_st2_into_st1_f]
        max_=maximum(dkappa_f_singleb4 .+ du_f_married_fv_singleb4_view .+ dEU_st1_f_view)
        num_p_st1_f_singleb4=exp.(dkappa_f_singleb4 .+ du_f_married_fv_singleb4_view .+ dEU_st1_f_view .- max_)
        den_p_st1_f_singleb4=sum(num_p_st1_f_singleb4,dims=(1,2,3)) .+ 1
        p_st1_f_singleb4=num_p_st1_f_singleb4./den_p_st1_f_singleb4
        # maximum(abs.(p_st1_f_singleb4.-CCP_true.p_st1_f_singleb4)./CCP_true.p_st1_f_singleb4)
        p0_st1_f_singleb4=  dropdims(1 ./den_p_st1_f_singleb4,dims=(1,2,3))
    end

    if i in [0, 4] 
        # Second stage structural CCPs
        du_m_single=sum(Xstruct_grid.Xstruct_st2_single .* bccp_st2_m_single,dims=1)
        du_m_single=reshape(du_m_single,(xbin_m,hbin_m,NTypes_m,J_m+1,T2))
        # To avoid overflow at low precisions
        max_=maximum(du_m_single   .+ fvt1_RX1.dEU_st2_single_m)
        num_p_st2_m_single=exp.(du_m_single   .+ fvt1_RX1.dEU_st2_single_m .-max_)
        den_p_st2_m_single=sum(num_p_st2_m_single,dims=4) 
        p_st2_m=num_p_st2_m_single./den_p_st2_m_single
        # maximum(abs.(p_st2_m.-CCP_true.p_st2_m)./CCP_true.p_st2_m)
        p0_st2_m=p_st2_m[:,:,:,1,:]

        du_f_single=sum(Xstruct_grid.Xstruct_st2_single .* bccp_st2_f_single,dims=1)
        du_f_single=reshape(du_f_single,(xbin_m,hbin_m,NTypes_m,J_m+1,T2))
        max_=maximum(du_f_single   .+ fvt1_RX1.dEU_st2_single_f)
        num_p_st2_f_single=exp.(du_f_single   .+ fvt1_RX1.dEU_st2_single_f .- max_)
        den_p_st2_f_single=sum(num_p_st2_f_single,dims=4) 
        p_st2_f=num_p_st2_f_single./den_p_st2_f_single
        # maximum(abs.(p_st2_f.-CCP_true.p_st2_f)./CCP_true.p_st2_f)
        p0_st2_f=p_st2_f[:,:,:,1,:]

        dv_m_married=du_m_married .+ fvt1_RX1.dEU_st2_married_m
        dv_f_married=du_f_married .+ fvt1_RX1.dEU_st2_married_f

        max_=maximum(Lambda .* dv_m_married .+ (1 .- Lambda) .* dv_f_married)
        num_p_st2_cpl= exp.(Lambda .* dv_m_married .+ (1 .- Lambda) .* dv_f_married .- max_)
        den_p_st2_cpl= sum(num_p_st2_cpl,dims=(7,8))
        p_st2_cpl=num_p_st2_cpl ./ den_p_st2_cpl
        # maximum(abs.(p_st2_cpl.-CCP_true.p_st2_cpl)./CCP_true.p_st2_cpl)
        p0_st2_cpl=p_st2_cpl[:,:,:,:,:,:,1,1,:]
    end


    CCP_struct= CCP_Structural(
        p0_st1_m_marriedb4, 
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
        p_st2_cpl
    )



    typeof(CCP_struct.p_st2_cpl)




    return CCP_struct
end


mutable struct CCP_Structural
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

