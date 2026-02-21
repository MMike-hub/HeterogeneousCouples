function resCCP(dt,CCP_,i=0)
    global dummy_size0

    res_st1_m_singleb4_data=dummy_size0 
    res_st1_f_singleb4_data=dummy_size0 
    res_st1_m_marriedb4_data=dummy_size0 
    res_st1_f_marriedb4_data=dummy_size0 
    res_st2_f_data=dummy_size0
    res_st2_m_data=dummy_size0 
    res_st2_cpl_data=dummy_size0
    
    if (i==1) | (i==0)
        p_st1_m_singleb4_data=@view CCP_.p_st1_m_singleb4[:,:,:, dt.partidx_p_st1_m_singleb4_data] 
        p_st1_m_singleb4_data=permutedims(p_st1_m_singleb4_data,(4,1,2,3))
        res_st1_m_singleb4_data=reshape(dt.I_st1_ .- p_st1_m_singleb4_data,(:,phibin_sp_dest*hbin_sp_dest*NTypes_sp_dest))
        
        p_st1_m_marriedb4_data=@view CCP_.p_st1_m_marriedb4[:,:,:,dt.partidx_p_st1_m_marriedb4_data] 
        p_st1_m_marriedb4_data=permutedims(p_st1_m_marriedb4_data,(4,1,2,3)) 
        res_st1_m_marriedb4_data=reshape(dt.I_st1_ .- p_st1_m_marriedb4_data,(:,phibin_sp_dest*hbin_sp_dest*NTypes_sp_dest))
    end
    if (i==2) | (i==0)
        p_st1_f_singleb4_data=@view CCP_.p_st1_f_singleb4[:,:,:,dt.partidx_p_st1_f_singleb4_data] 
        p_st1_f_singleb4_data=permutedims(p_st1_f_singleb4_data,(4,1,2,3))
        res_st1_f_singleb4_data=reshape(dt.I_st1_ .- p_st1_f_singleb4_data,(:,phibin_sp_dest*hbin_sp_dest*NTypes_sp_dest))

        p_st1_f_marriedb4_data=@view CCP_.p_st1_f_marriedb4[:,:,:,dt.partidx_p_st1_f_marriedb4_data] 
        p_st1_f_marriedb4_data=permutedims(p_st1_f_marriedb4_data,(4,1,2,3)) 
        res_st1_f_marriedb4_data=reshape(dt.I_st1_ .- p_st1_f_marriedb4_data,(:,phibin_sp_dest*hbin_sp_dest*NTypes_sp_dest))
    end
    if (i==3) | (i==0)
        p_st2_m_data=@view permutedims(CCP_.p_st2_m,(4,1,2,3,5))[:,dt.partidx_p_st2_m_data] 
        p_st2_m_data= permutedims(p_st2_m_data,(2,1))
        res_st2_m_data=(dt.I_st2_single_ .- p_st2_m_data)[:,2:end]
    end
    if (i==4) | (i==0)
        p_st2_f_data=@view permutedims(CCP_.p_st2_f,(4,1,2,3,5))[:,dt.partidx_p_st2_f_data] 
        p_st2_f_data=permutedims(p_st2_f_data,(2,1))
        res_st2_f_data=(dt.I_st2_single_ .- p_st2_f_data)[:,2:end] #Drop the first residual, corresponding to the baseline choice, otherwise residuals are linearly dependent
    end
    if (i==5) | (i==0)
        p_st2_cpl_data=@view permutedims(CCP_.p_st2_cpl,(7,8,1,2,3,4,5,6,9))[:,:,dt.partidx_p_st2_cpl_data]
        p_st2_cpl_data=permutedims(p_st2_cpl_data,(3,1,2))
        res_st2_cpl_data=reshape(dt.I_st2_cpl_ .- p_st2_cpl_data,(:,(J_m+1)*(J_f+1)))[:,2:end]
    end 

    res=Residuals_(res_st1_m_singleb4_data, res_st1_f_singleb4_data, res_st1_m_marriedb4_data, res_st1_f_marriedb4_data, res_st2_f_data, res_st2_m_data, res_st2_cpl_data)

    return res
end

struct Residuals_
    res_st1_m_singleb4_data
    res_st1_f_singleb4_data
    res_st1_m_marriedb4_data
    res_st1_f_marriedb4_data
    res_st2_f_data
    res_st2_m_data
    res_st2_cpl_data
end