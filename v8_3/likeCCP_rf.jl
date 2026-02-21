# This speficiation of the likelihood allows the calculation of the likelihood separately for men, women, 1st stage and 2nd stage
function likeCCP_rf(dt,CCP_,i=0)
    global dummy_sizeNT, dummy_size0
    # Initialize CCPs as views (no new memory allocation) of a vector of ones
    dt_p_st1_m_singleb4=  dummy_sizeNT
    dt_p0_st1_m_singleb4=  dummy_sizeNT
    dt_p_st1_m_marriedb4=  dummy_sizeNT
    dt_p0_st1_m_marriedb4=  dummy_sizeNT
    dt_p_st1_f_singleb4=  dummy_sizeNT
    dt_p0_st1_f_singleb4=  dummy_sizeNT
    dt_p_st1_f_marriedb4=  dummy_sizeNT
    dt_p0_st1_f_marriedb4=  dummy_sizeNT
    dt_p_st2_m_single=  dummy_sizeNT
    dt_p_st2_f_single=  dummy_sizeNT
    dt_p_st2_cpl=  dummy_sizeNT


    if (i==1) | ( i==0)
        # Because in the regression errors I have one observation for each combination of unobserved types of agent and spouse, and because single agents don't have a spouse, CCPs corresponding to first-stage choices of people who choose to be single appear twice identically in the Likelihood. For this reason, I take their square root.
        if (i==1)
            # Select only the relevant observations
            if length(dt.idx_p_st1_m_singleb4_reduced)>1
                dt_p_st1_m_singleb4= @view CCP_.p_st1_m_singleb4[dt.idx_p_st1_m_singleb4_reduced]
            else
                dt_p_st1_m_singleb4=dummy_size0
            end
            if length(dt.idx_p0_st1_m_singleb4_reduced)>1
                dt_p0_st1_m_singleb4= @view CCP_.p0_st1_m_singleb4[dt.idx_p0_st1_m_singleb4_reduced]
                dt_p0_st1_m_singleb4=sqrt.(dt_p0_st1_m_singleb4)
            else
                dt_p0_st1_m_singleb4=dummy_size0
            end
            if length(dt.idx_p_st1_m_marriedb4_reduced)>1
                dt_p_st1_m_marriedb4= @view CCP_.p_st1_m_marriedb4[dt.idx_p_st1_m_marriedb4_reduced]
            else
                dt_p_st1_m_marriedb4=dummy_size0
            end
            if length(dt.idx_p0_st1_m_marriedb4_reduced)>1
                dt_p0_st1_m_marriedb4= @view CCP_.p0_st1_m_marriedb4[dt.idx_p0_st1_m_marriedb4_reduced]
                dt_p0_st1_m_marriedb4=sqrt.(dt_p0_st1_m_marriedb4)
            else
                dt_p0_st1_m_marriedb4=dummy_size0
            end
            Like=vcat(dt_p_st1_m_singleb4,dt_p0_st1_m_singleb4,dt_p_st1_m_marriedb4,dt_p0_st1_m_marriedb4)
        else
            dt_p_st1_m_singleb4= @view CCP_.p_st1_m_singleb4[dt.idx_p_st1_m_singleb4]
            dt_p0_st1_m_singleb4= @view CCP_.p0_st1_m_singleb4[dt.idx_p0_st1_m_singleb4]
            dt_p_st1_m_marriedb4= @view CCP_.p_st1_m_marriedb4[dt.idx_p_st1_m_marriedb4]
            dt_p0_st1_m_marriedb4= @view CCP_.p0_st1_m_marriedb4[dt.idx_p0_st1_m_marriedb4]

            Like_st1_marry_m=(dt_p_st1_m_singleb4.*(dt.married_tminus1_.==0) .+ dt_p_st1_m_marriedb4 .* dt.married_tminus1_ ) .*(dt.sex_.==1)
            Like_st1_single_m=(dt_p0_st1_m_singleb4.*(dt.married_tminus1_.==0) .+ dt_p0_st1_m_marriedb4 .* dt.married_tminus1_ ) .*(dt.sex_.==1)
            Like_st1_single_m=sqrt.(Like_st1_single_m)
            Like_st1_m=Like_st1_marry_m .* dt.married_ .+ (dt.married_.==0) .* Like_st1_single_m
            # The first observation of a household does not contain first-step probabilities
            # It's a contrived way to replace Like_st1 with ones when first sampled, but Like_st1[dt.firstampled.==0].=1 fails with Zygote
            Like_st1_m= dummy_sizeNT .* dt.firstsampled .+ Like_st1_m .* (dt.firstsampled.==0)
        end
    end
    if (i==2) | ( i==0)
        if (i==2)
            # Select only the relevant observations
            if length(dt.idx_p_st1_f_singleb4_reduced)>1
                dt_p_st1_f_singleb4= @view CCP_.p_st1_f_singleb4[dt.idx_p_st1_f_singleb4_reduced]
            else
                dt_p_st1_f_singleb4=dummy_size0
            end
            if length(dt.idx_p0_st1_f_singleb4_reduced)>1
                dt_p0_st1_f_singleb4= @view CCP_.p0_st1_f_singleb4[dt.idx_p0_st1_f_singleb4_reduced]
                dt_p0_st1_f_singleb4=sqrt.(dt_p0_st1_f_singleb4)
            else
                dt_p0_st1_f_singleb4=dummy_size0
            end
            if length(dt.idx_p_st1_f_marriedb4_reduced)>1
                dt_p_st1_f_marriedb4= @view CCP_.p_st1_f_marriedb4[dt.idx_p_st1_f_marriedb4_reduced]
            else
                dt_p_st1_f_marriedb4=dummy_size0
            end
            if length(dt.idx_p0_st1_f_marriedb4_reduced)>1
                dt_p0_st1_f_marriedb4= @view CCP_.p0_st1_f_marriedb4[dt.idx_p0_st1_f_marriedb4_reduced]
                dt_p0_st1_f_marriedb4=sqrt.(dt_p0_st1_f_marriedb4)
            else
                dt_p0_st1_f_marriedb4=dummy_size0
            end
            Like=vcat(dt_p_st1_f_singleb4,dt_p_st1_f_marriedb4,dt_p0_st1_f_singleb4,dt_p0_st1_f_marriedb4)
        else
            dt_p_st1_f_singleb4= @view CCP_.p_st1_f_singleb4[dt.idx_p_st1_f_singleb4]
            dt_p0_st1_f_singleb4= @view CCP_.p0_st1_f_singleb4[dt.idx_p0_st1_f_singleb4]
            dt_p_st1_f_marriedb4= @view CCP_.p_st1_f_marriedb4[dt.idx_p_st1_f_marriedb4]
            dt_p0_st1_f_marriedb4= @view CCP_.p0_st1_f_marriedb4[dt.idx_p0_st1_f_marriedb4]

            Like_st1_marry_f=(dt_p_st1_f_singleb4.*(dt.married_tminus1_.==0) .+ dt_p_st1_f_marriedb4 .* dt.married_tminus1_ ) .*(dt.sex_.==2)
            Like_st1_single_f=(dt_p0_st1_f_singleb4.*(dt.married_tminus1_.==0) .+ dt_p0_st1_f_marriedb4 .* dt.married_tminus1_ ) .*(dt.sex_.==2)
            Like_st1_single_f=sqrt.(Like_st1_single_f)
            Like_st1_f=Like_st1_marry_f .* dt.married_ .+ (dt.married_.==0) .* Like_st1_single_f
            Like_st1_f= dummy_sizeNT .* dt.firstsampled .+ Like_st1_f .* (dt.firstsampled.==0)
        end
    end
    if (i==3) | ( i==0)
        # Because in the regression errors I have one observation for each combination of unobserved types of agent and spouse, and because single agents don't have a spouse, CCPs corresponding to second-stage choices of single people appear twice in the Likelihood. For this reason, I take their square root.
        if (i==3)
            if length(dt.idx_p_st2_m_single_reduced)>1
                Like= @view CCP_.p_st2_m[dt.idx_p_st2_m_single_reduced]
                Like=sqrt.(Like)
            else
                Like=0
            end
        else
            dt_p_st2_m_single= @view CCP_.p_st2_m[dt.idx_p_st2_m_single]
            Like_st2_m=dt_p_st2_m_single .* (dt.sex_.==1) .* (dt.married_.==0)
            Like_st2_m=sqrt.(Like_st2_m)
        end
    end
    if (i==4) | ( i==0)
        if (i==4)
            if length(dt.idx_p_st2_f_single_reduced)>1
                Like= @view CCP_.p_st2_f[dt.idx_p_st2_f_single_reduced]
                Like=sqrt.(Like)
            else
                Like=0
            end
        else
            dt_p_st2_f_single= @view CCP_.p_st2_f[dt.idx_p_st2_f_single]
            Like_st2_f=dt_p_st2_f_single .* (dt.sex_.==2) .* (dt.married_.==0)
            Like_st2_f=sqrt.(Like_st2_f)
        end
    end
    if (i==5) | ( i==0)
        # To avoid entering twice the second-stage choice of the originally sampled couples, I should take the square root of Like_st2
        if (i==5)
            if length(dt.idx_p_st2_cpl_reduced)>1
                Like= @view CCP_.p_st2_cpl[dt.idx_p_st2_cpl_reduced]
                if sum(dt_rf.originalspouse_reduced)>0
                    Like=sqrt.(Like).*dt.originalspouse_reduced .+ (1 .- dt.originalspouse_reduced).*Like
                end
            else
                Like=0
            end
        else
            dt_p_st2_cpl= @view CCP_.p_st2_cpl[dt.idx_p_st2_cpl]
            Like_st2_cpl=dt_p_st2_cpl .* dt.married_
            if sum(dt.originalspouse)>0
                Like_st2_cpl= sqrt.(Like_st2_cpl).*dt.originalspouse .+ (1 .- dt.originalspouse).*Like_st2_cpl
            end
        end
    end

    
    if (i==0)
        Like =  (Like_st1_m .* (dt.sex_.==1) .+ Like_st1_f .* (dt.sex_.==2)) .* (Like_st2_m .+ Like_st2_f .+ Like_st2_cpl)

        # Initial p(x,z|s):
        # pinit_data=intcond_nonparametric(dt_rf,dt,PType)
        # pinit_data=ones(precision, size(Like))
        # pinit_data=dummy
        # Like=Like.*pinit_data

        # Like = ifelse.(dt.insample.==0,missing,Like)
        # Like[dt.insample!_Bool] .= 1
        # This more contrived way is compatible with Zygote
        Like=(1 .- dt.insample) .+ dt.insample .* Like
    end
    return Like
end

function intcond_nonparametric(dt_rf,dt,PType)
    I_intcond_single=zeros(precision,size(dt.phi,1),phibin,hbin)
    married_=coalesce.(dt.married,0)
    for phi in 1:phibin
        for h in 1:hbin
            I_intcond_single[:,phi,h]=prod(dt.phi.===phival[phi],dims=2).*(dt.h.===hval[h]).*(dt.firstsampled.==1).*(married_.==0)
        end
    end
    I_intcond_single=reshape(I_intcond_single.*PType,(2*N*T2,NTypes^2,phibin,hbin))
    denominator=reshape((dt.firstsampled.==1).*(married_.==0),(2*N*T2,NTypes^2))
    pinit_single=sum(I_intcond_single,dims=1)./sum(denominator, dims=1)
    pinit_single=dropdims(pinit_single,dims=1)

    I_intcond_married=zeros(precision,size(dt.phi,1),phibin,phibin_sp,hbin,hbin_sp)
    for phi in 1:phibin
        for h in 1:hbin
            for phi_sp in 1:phibin_sp
                for h_sp in 1:hbin_sp
                    I_intcond_married[:,phi,phi_sp,h,h_sp]=prod(dt.phi.===phival[phi],dims=2).*(dt.h.===hval[h]).*
                                                        prod(dt.phi_sp_dest.===phival[phi_sp],dims=2).*(dt.h_sp_dest.===hval[h_sp]).*(dt.firstsampled.==1).*(married_.==1)
                end
            end
        end
    end
    I_intcond_married=reshape(I_intcond_married.*PType,(2*N*T2,NTypes^2,phibin,phibin_sp,hbin,hbin_sp))
    denominator=reshape((dt.firstsampled.==1).*(married_.==1),(2*N*T2,NTypes^2))
    pinit_married=sum(I_intcond_married,dims=1)./sum(denominator, dims=1)
    pinit_married=dropdims(pinit_married,dims=1)

    idx_single=CartesianIndex.(coalesce.(dt.s,1),coalesce.(dt.phistate,1),coalesce.(dt.hstate,1))
    idx_married=CartesianIndex.(coalesce.(dt.s,1),coalesce.(dt.phistate,1),coalesce.(dt.phistate_sp_dest,1),coalesce.(dt.hstate,1),coalesce.(dt.hstate_sp_dest,1))
    pinit_data=ones(precision, size(dt.x)).*(1 .-dt.firstsampled) .+ (pinit_single[idx_single].*(married_==0) .+ pinit_married[idx_married].*(married_==1)).*dt.firstsampled

    return pinit_data
end