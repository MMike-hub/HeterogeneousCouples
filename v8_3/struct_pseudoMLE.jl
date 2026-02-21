# THE STRUCTURAL LOGIT IS CHARACTERIZED BY PARAMETER VALUES THAT ARE MOSTLY THE SAME ACROSS CHOICES, WHILE IT'S THE RELEVANT REGRESSORS
# THAT CHANGE ACROSS CHOICES. No reshape of the parameter vector needed here.
# Unlike the reduced form logit, where all parameter values necessarily change across choices 
# and consequently a reshape of the parameter vector is needed

function struct_pseudoMLE_closure(dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i=0, bccp_dict_pseudoMLE=[])
    return bccp_ -> f_struct_pseudoMLE(bccp_, dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i, bccp_dict_pseudoMLE)
end

function f_struct_pseudoMLE(bccp_,dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i=0, bccp_dict_pseudoMLE=[])
    CCP_struct=CCP_structural(bccp_ ,Lambda, Xstruct_grid, fvt1_RX1, CCP_, i, bccp_dict_pseudoMLE)
    Like = likeCCP_struct(dt, CCP_struct, i)

    if i in [1,3]
        if length(dt.idx_p_st1_m_singleb4_reduced_obs)>1
            PType_p_st1_m_singleb4= @view PType[dt.idx_p_st1_m_singleb4_reduced_obs]
        else
            PType_p_st1_m_singleb4=dummy_size0
        end
        if length(dt.idx_p0_st1_m_singleb4_reduced_obs)>1
            PType_p0_st1_m_singleb4= @view PType[dt.idx_p0_st1_m_singleb4_reduced_obs]
        else
            PType_p0_st1_m_singleb4=dummy_size0
        end
        if length(dt.idx_p_st1_m_marriedb4_reduced_obs)>1
            PType_p_st1_m_marriedb4= @view PType[dt.idx_p_st1_m_marriedb4_reduced_obs]
        else
            PType_p_st1_m_marriedb4=dummy_size0
        end
        if length(dt.idx_p0_st1_m_marriedb4_reduced_obs)>1
            PType_p0_st1_m_marriedb4= @view PType[dt.idx_p0_st1_m_marriedb4_reduced_obs]
        else
            PType_p0_st1_m_marriedb4=dummy_size0
        end
        PType_st1_m=vcat(PType_p_st1_m_singleb4, PType_p0_st1_m_singleb4, PType_p_st1_m_marriedb4,PType_p0_st1_m_marriedb4)
        if (i==1)
            PType_=PType_st1_m
        end
    end

    if i in [2,3]
        if length(dt.idx_p_st1_f_singleb4_reduced_obs)>1
            PType_p_st1_f_singleb4= @view PType[dt.idx_p_st1_f_singleb4_reduced_obs]
        else
            PType_p_st1_f_singleb4=dummy_size0
        end
        if length(dt.idx_p0_st1_f_singleb4_reduced_obs)>1
            PType_p0_st1_f_singleb4= @view PType[dt.idx_p0_st1_f_singleb4_reduced_obs]
        else
            PType_p0_st1_f_singleb4=dummy_size0
        end
        if length(dt.idx_p_st1_f_marriedb4_reduced_obs)>1
            PType_p_st1_f_marriedb4= @view PType[dt.idx_p_st1_f_marriedb4_reduced_obs]
        else
            PType_p_st1_f_marriedb4=dummy_size0
        end
        if length(dt.idx_p0_st1_f_marriedb4_reduced_obs)>1
            PType_p0_st1_f_marriedb4= @view PType[dt.idx_p0_st1_f_marriedb4_reduced_obs]
        else
            PType_p0_st1_f_marriedb4=dummy_size0
        end
        PType_st1_f=vcat(PType_p_st1_f_singleb4, PType_p0_st1_f_singleb4, PType_p_st1_f_marriedb4,PType_p0_st1_f_marriedb4)
        if (i==2)
            PType_=PType_st1_f
        end
        if (i==3)
            PType_=vcat(PType_st1_m,PType_st1_f)
        end
    end

    if (i==4)
        if length(dt.idx_p_st2_m_single_reduced_obs)>1
            PType_st2_m_single= @view PType[dt.idx_p_st2_m_single_reduced_obs]
        else
            PType_st2_m_single=dummy_size0
        end
        if length(dt.idx_p_st2_f_single_reduced_obs)>1
            PType_st2_f_single= @view PType[dt.idx_p_st2_f_single_reduced_obs]
        else
            PType_st2_f_single=dummy_size0
        end
        if length(dt.idx_p_st2_cpl_reduced_obs)>1
            PType_st2_cpl= @view PType[dt.idx_p_st2_cpl_reduced_obs]
        else
            PType_st2_cpl=dummy_size0
        end
        PType_=vcat(PType_st2_m_single,PType_st2_f_single,PType_st2_cpl)
    end

    if (i==0)
        PType_=PType
    end
    
    Like = PType_' * log.(Like)
    # This is just to ensure that Like is a scalar
    Like=-sum(Like)
    return Like
end

