function rf_pseudoMLE_closure(dt, dt_rf_grid, dt_rf, PType, i=0)
    # If using analytic derivatives
    # return (F,G,b1_) -> wlogit_rf_fg!(F,G,b1_, dt, dt_rf_grid, dt_rf, PType, i)

    # If using Optim's automatic differentiation or providing a gradient to optim separately
    return b1_ ->wlogit_rf(b1_, dt, dt_rf_grid, dt_rf, PType, i)
end

function wlogit_rf_fg!(F,G,b1_, dt, dt_rf_grid, dt_rf, PType, i=0)
    global b1_previous_inner
    diff=abs(sum(b1_.-b1_previous_inner))
    println("b1_change=$diff")
    b1_previous_inner=b1_

    CCP_rf=CCP_reducedform(b1_, dt_rf_grid, dt_rf, i)
    println("CCP_rf updated")
    if G !== nothing
        G.=-likeCCP_rf_gradient(dt,CCP_rf,PType,i)
        println("gradient updated")
    end
    if F !== nothing
        Like=likeCCP(dt,CCP_rf,i)
        Like=-sum(PType' * log.(Like))
        println("objective updated")
        return Like
    end
end

function wlogit_rf(b1_, dt, dt_rf_grid, dt_rf, PType, i=0)
    global dummy_size0
    CCP_rf=CCP_reducedform(b1_, dt_rf_grid, dt_rf, i)
    Like=likeCCP_rf(dt, CCP_rf, i)

    if (i==1)
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
        PType_=vcat(PType_p_st1_m_singleb4, PType_p0_st1_m_singleb4, PType_p_st1_m_marriedb4, PType_p0_st1_m_marriedb4)
    elseif (i==2)
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
        PType_=vcat(PType_p_st1_f_singleb4, PType_p0_st1_f_singleb4, PType_p_st1_f_marriedb4,PType_p0_st1_f_marriedb4)
    elseif (i==3)
        if length(dt.idx_p_st2_m_single_reduced_obs)>1
            PType_= @view PType[dt.idx_p_st2_m_single_reduced_obs]
        else
            PType_=0
        end
    elseif (i==4)
        if length(dt.idx_p_st2_f_single_reduced_obs)>1
            PType_= @view PType[dt.idx_p_st2_f_single_reduced_obs]
        else
            PType_=0
        end
    elseif (i==5)
        if length(dt.idx_p_st2_cpl_reduced_obs)>1
            PType_= @view PType[dt.idx_p_st2_cpl_reduced_obs]
        else
            PType_=0
        end
    elseif (i==0)
        PType_=PType
    end

    Like=-sum(PType_' * log.(Like))
    return Like
end