function likeCCP_rf_gradient(dt,CCP_,PType,i=0)
    # Initialize CCPs as views (no new memory allocation) of a vector of ones

    global dummy_size0, dummy_size10

    d_Like_st1_m=dummy_size10
    d_Like_st1_f=dummy_size10
    d_Like_st2_m=dummy_size10
    d_Like_st2_f=dummy_size10
    d_Like_st2_cpl=dummy_size10

    if (i==1) | ( i==0)
        view1= @view CCP_.p_st1_m_singleb4[dt_rf.idx_st1_m_singleb4_allchoices_observedstate]
        dt_p_st1_m_singleb4=  ( (dt_rf.idx_st1_m_allchoices.==dt.idx_p_st1_m_singleb4_choicesonly).* (dt.married_.==1)  .- view1)
        view2= @view CCP_.p_st1_m_marriedb4[dt_rf.idx_st1_m_marriedb4_allchoices_observedstate]
        dt_p_st1_m_marriedb4= ( (dt_rf.idx_st1_m_allchoices.==dt.idx_p_st1_m_marriedb4_choicesonly).* (dt.married_.==1)  .- view2)

        d_Like_st1_m=reshape((dt_p_st1_m_singleb4.*(dt.married_tminus1_.==0) .+ dt_p_st1_m_marriedb4.* dt.married_tminus1_) ,(:,dt_rf.dim_x_st1_2,dt_rf.dim_y_st1)) .* dt_rf.xx_st1 .* (dt.sex_.==1)  .* dt.insample .* PType
        d_Like_st1_m=sum(reshape(d_Like_st1_m,(:,dt_rf.dim_y_st1*dt_rf.dim_x_st1_2)), dims=1)
    end
    if (i==2) | ( i==0)
        view1= @view CCP_.p_st1_f_singleb4[dt_rf.idx_st1_f_singleb4_allchoices_observedstate]
        dt_p_st1_f_singleb4=  ( (dt_rf.idx_st1_f_allchoices.==dt.idx_p_st1_f_singleb4_choicesonly).* (dt.married_.==1)  .- view1)
        view2= @view CCP_.p_st1_f_marriedb4[dt_rf.idx_st1_f_marriedb4_allchoices_observedstate]
        dt_p_st1_f_marriedb4= ( (dt_rf.idx_st1_f_allchoices.==dt.idx_p_st1_f_marriedb4_choicesonly).* (dt.married_.==1)  .- view2)

        d_Like_st1_f=reshape((dt_p_st1_f_singleb4.*(dt.married_tminus1_.==0) .+ dt_p_st1_f_marriedb4.* dt.married_tminus1_) ,(:,dt_rf.dim_x_st1_2,dt_rf.dim_y_st1)) .* dt_rf.xx_st1 .* (dt.sex_.==1)  .* dt.insample .* PType
        d_Like_st1_f=sum(reshape(d_Like_st1_f,(:,dt_rf.dim_y_st1*dt_rf.dim_x_st1_2)), dims=1)
    end
    if (i==3) | ( i==0)
        dt_p_st2_m_single= @view CCP_.p_st2_m[dt_rf.idx_st2_f_allchoices_observedstate]
    end
    if (i==4) | ( i==0)
        dt_p_st2_f_single= @view CCP_.p_st2_f[dt.idx_p_st2_f_single]
    end
    if (i==5) | ( i==0)
        dt_p_st2_cpl= @view CCP_.p_st2_cpl[dt.idx_p_st2_cpl]
    end

    
    # Like_st1_f_marry=(dt_p_st1_f_singleb4.*(1 .- dt.married_tminus1) .+ dt_p_st1_f_marriedb4 .* dt.married_tminus1 ) .*(dt.sex.===2)
    # Like_st1_f_single=(dt_p0_st1_f_singleb4.*(1 .- dt.married_tminus1) .+ dt_p0_st1_f_marriedb4 .* dt.married_tminus1 ) .*(dt.sex.===2)
    
    
    # Like_st1=Like_st1_marry .*dt.married .+ (1 .- dt.married) .* Like_st1_single
    # Like_st1[Bool.(dt.firstsampled),:].=0 #The first observation of a household does not contain first-step probabilities
    
    # Like_st2= (dt_p_st2_m_single .* (1 .- dt.married) .* (dt.sex.===1) .+ dt_p_st2_f_single .* (1 .- dt.married) .* (dt.sex.===2) .+
    # dt_p_st2_cpl .* dt.married)

    # To avoid entering twice the second-stage choice of the originally sampled couples, I take the square root of Like_st2
     
    
    d_Like =  vcat(d_Like_st1_m[:], d_Like_st1_f[:], d_Like_st2_m[:], d_Like_st2_f[:], d_Like_st2_cpl[:])
        
    # Initial p(x,z|s):
    # pinit_data=intcond_nonparametric(dt_rf,dt,PType)
    # pinit_data=ones(precision,size(Like))
    # Like=Like.*pinit_data

    # Sanity check.
    # Like2=coalesce.(reshape(Like,(2*N,T2,NTypes^2)),0)
    # NotMarried=coalesce.(1 .-vcat(dt.Married_1,dt.Married_2),0)
    # check1= maximum(Like2[:,:,1].*NotMarried .- Like2[:,:,3].*NotMarried) > eps()*10
    # check2 = maximum(Like2[:,:,2].*NotMarried .- Like2[:,:,4].*NotMarried) > eps()*10
    # if check1==true | check2==true
    #     println("ISSUE HERE: likeCCP")
    # end

    return d_Like
end



