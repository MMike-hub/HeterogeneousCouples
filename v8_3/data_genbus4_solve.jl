# Bus data generation
    # Backwards induction to solve model, i.e., find the continuation values big FV and FU
    # Notice that the lifespan is assumed to be T (column T+1 in FV is all zeros, i.e., utility after last period is zero)

    function f_data_genbus4_solve()
        """
        Generates decision data for a CCP Model
            alpha [vector float]: Vector of coefficients to parameterize the model 
            N_ [int]: Number of individuals 
            T [int]: Number of periods to simulate 
            xtran [Matrix float]: Transition probability matrix
            xval [int]: number of observed states
        
        Returns: 
            Y [matrix int]: NxT matrix of decisions of each individual at time t  
            Phi [matrix float]: NxT matrix of observed states of each individual at time t
            Phistate [matrix int]: NxT matrix of index of states of each individual at time t
            State [matrix int]: NxT matrix of unobs. state of each individual at time t
            FV [matrix float]: FOLLOWING OUR PAPER'S NOTATION. FV is the expected value at stage 2,  i.e. after the marriage decision is taken and before the other choice shocks are realized.
            FU [array float]: FOLLOWING OUR PAPER'S NOTATION. FU is the expected value at stage 1, i.e. before the marriage decision
        """
    
        # alpha=parameters
        # N_=sample size
        # T=time horizon
        # xtran=transition matrix
        # xtranc=cumulative transition probability, also reshaped with respect to xtran
        # xbin, hbin=size of state space in each dimension
        # xval, hval= discretized state space dimensions
        Beta = alpha[5]
        # Pi_m=proportion of males of unobserved type s=1
        Pi_m = alpha[6]
        Pi_f = alpha[7]
        eul = 0.577215665

    
        ##########################
        # Second-stage objects
        ##########################
        FV = zeros(xbin,hbin, NTypes, NSex, T  + 1)
        FV_married=zeros(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_f, NSex, T  + 1)
        vJ=zeros(xbin,hbin,NTypes,J+1,NSex,T)
        p_st2=zeros(xbin,hbin,NTypes,J+1,NSex,T)
        vJ_married=zeros(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_f,J_m+1,J_f+1,NSex,T)
        # p_st2_cpl gives the probability that a couple of type (phi,x_sp,h,h_sp,s,s_sp) at time t takes a discrete decision (j,j_sp)
        p_st2_cpl=zeros(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_f,J_m+1,J_f+1,T)
        lambda=permutedims(Lambda_T,(1,4,2,5,3,6,7,8,9)) # rearrange dimensions as (xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_f,1,1,T) # The 2 placeholder dimensions correspond to J_m+1 and J_f+1
        lambda=reshape(lambda,(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_f,1,1,1,T)) # Add one placeolder dimension (sex), one for male action j_m and one for female action j_f for conformability to element-wise multiplication with vJ
    
        
        ##########################
        # First-stage objects
        ########################## 
    
        # Due to transition costs in and out of marriage, the expected value at the beginning of stage 1 also depends on the state of the partner.
        # FU contains the expected value of single agents at the beginning of stage 1 
        FU = zeros(phibin, hbin, NTypes, NSex, T  + 1)
        # FU_married contains the expected value of married agents at the beginning of stage 1, as a function of who they are married to
        FU_married=zeros(phibin,phibin_sp_orig,hbin, hbin_sp_orig, NTypes, NTypes_sp_orig, NSex, T  + 1)
        # u_ is the conditional valuation function of the marriage choices at stage 1
        u_marry_marriedb4=zeros(phibin,phibin_sp_orig,phibin_sp_dest,hbin,hbin_sp_orig,hbin_sp_dest,NTypes,NTypes_sp_orig,NTypes_sp_dest,NSex,T) # The dimensions represent (phi,phi_spouse_origin,phi_spouse_destination,h,h_spouse_origin,h_spouse_destination,s,s_spouse_origin,s_spouse_destination,sex)
        p_st1_marriedb4=zeros(phibin,phibin_sp_orig,phibin_sp_dest,hbin,hbin_sp_orig,hbin_sp_dest,NTypes,NTypes_sp_orig,NTypes_sp_dest,NSex,T)
        u_single_marriedb4=zeros(phibin,phibin_sp_orig,hbin,hbin_sp_orig,NTypes,NTypes_sp_orig,NSex,T) # The dimensions represent (phi,phi_spouse_origin,h,h_spouse_origin,s,s_spouse_origin,sex)
        p0_st1_marriedb4=zeros(phibin,phibin_sp_orig,hbin,hbin_sp_orig,NTypes,NTypes_sp_orig,NSex,T)
        u_marry_single=zeros(phibin,phibin_sp_dest,hbin,hbin_sp_dest,NTypes,NTypes_sp_dest,NSex,T) # The dimensions represent (phi,phi_spouse_destination,h,h_spouse_destination,s,s_spouse_destination,sex)
        p_st1_singleb4=zeros(phibin,phibin_sp_dest,hbin,hbin_sp_dest,NTypes,NTypes_sp_dest,NSex,T)
        u_single_single=zeros(phibin,hbin,NTypes,NSex,T) # The dimensions represent (phi,h,s,sex)
        p0_st1_singleb4=zeros(phibin,hbin,NTypes,NSex,T) 


        # The kappa function computes transition costs.
        function kappa(phi,phi_sp_orig,phi_sp_dest,h,h_sp_orig,h_sp_dest,s,s_sp_orig,s_sp_dest,sex)
            # Establish normalization of 1st stage flow utilities
            # If COND_ORIG_SPOUSE==false, transition costs cannot depend on characteristics of the previous spouse
            # Instead, transition cost to being single should always be equal to 0 or the model is not identified.
            local kappa
            if phi_sp_dest>0
                kappa=  bk[1] #+ bk[1] .* phival[phi_sp_dest,2]  
            else
                kappa=0
            end
            return kappa
        end
    
        function util_fn(x,h,s,j,sex,t)
            x_=xval[x]
            h_=hval[h]
            if j==1
                # Establish normalization of 2nd stage flow utilities
                util=0
            else
                util=alpha[1]  + alpha[2] * x_ + alpha[3] * (s==2) + alpha[4] * j + Adj[t]
            end
            return util
        end
    
        function util_married_fn(x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,sex,t)
            local x_m_, h_m_, x_f_, h_f_
            global hval
            x_m_=xval[x_m]
            h_m_=hval[h_m]
            x_f_=xval[x_f]
            h_f_=hval[h_f]
            if sex==1 
                x_=x_m_
                s=s_m
                j=j_m
            else
                x_=x_f_
                s=s_f
                j=j_f
            end
            # Establish normalization of 2nd stage flow utilities
            if j==1
                util=0
            else
                util= alpha[1]  + alpha[2] * x_ + alpha[3] * (s==2) + alpha[4] * j + Adj[t]
            end
            return util
        end
    
    
        tic=time()
        lk = ReentrantLock()
        t = T 
        while t >= 1
            println(t)
            ##################################
            # Second-stage choices as single
            ##################################
            for s in 1:NTypes
                for h in 1:hbin
                    Threads.@threads for x in 1:xbin 
                        for sex=1:2
                            for j=1:J+1
                                j_tran=max(j-1,1) # Ensure that j=1 and j=2 both use xtran[x, :,h, 1]  
                                lock(lk) do 
                                    vJ[x,h,s,j,sex,t] = util_fn(x,h,s,j,sex, t) + xtran[x, :,h, j_tran]' * FU[:, h, s, sex, t+1]
                                end
                            end
                            lock(lk) do
                                # FV[x,h, s, sex, t] =  (log(sum(exp.(vJ[x,h,s,:,sex,t]))) + eul) 
                                FV[x,h, s, sex, t] =  log(sum(exp.(vJ[x,h,s,:,sex,t]))) 
                            end
                            for j=1:J+1
                                lock(lk) do
                                    p_st2[x,h,s,j,sex,t] = exp(vJ[x,h,s,j,sex,t])/(sum(exp.(vJ[x,h,s,:,sex,t])))
                                end
                            end
                            check=sum(p_st2[x,h,s,:,sex,t])
                            if abs(check-1)>1e-10
                                println("CCPs don't sum up to 1! $check")
                            end
                        end
                    end
                end
            end
            ######################################
            #  Second-stage choices as married
            ######################################
            for sex in 1:2
                for s_m in 1:NTypes_m
                    for s_f in 1:NTypes_f
                        for h_m in 1:hbin_m
                            for h_f in 1:hbin_f
                                if sex==1
                                    h=h_m
                                    h_sp=h_f
                                    s=s_m
                                    s_sp=s_f
                                else
                                    h=h_f
                                    h_sp=h_m
                                    s=s_f
                                    s_sp=s_m
                                end
                                if !COND_ORIG_SPOUSE
                                    h_sp=1
                                    s_sp=1
                                end
                                for x_m in 1:xbin_m
                                    # Threads.@threads for x_f in 1:xbin_f
                                    for x_f in 1:xbin
                                        for j_m=1:J+1
                                            for j_f=1:J+1
                                                lock(lk) do
                                                    # println("x_m is $x_m, x_f is $x_f, h_m is $h_m, h_f is $h_f, s_m is $s_m, s_f is $s_f, j_m is $j_f, sex is $sex, t is $t")
                                                    vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,sex,t] = util_married_fn(x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,sex,t) + sum(sum(xtran_married[:,:, x_m, x_f, h_m, h_f, j_m,j_f] .* FU_married[:,:, h, h_sp, s, s_sp, sex, t+1],dims=(1,2)))
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            for s_m in 1:NTypes_m
                for s_f in 1:NTypes_f
                    for h_m in 1:hbin_m
                        for h_f in 1:hbin_f
                            for x_m in 1:xbin_m
                                Threads.@threads for x_f in 1:xbin_f
                                # for x_f in 1:xbin_f
                                    denom=sum(exp.(lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t].*vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,1,t] .+ (1 .-lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t]).*vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,2,t]))
                                    num=exp.(lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t]* vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,1,t] + (1 - lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t]) * vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,2,t])
                                    lock(lk) do
                                        p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t].=num./denom
                                    end
                                    check=sum(p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t])
                                    if abs(check-1)>1e-15
                                        println("CCPs don't sum up to 1! $check")
                                    end
                                    denom=0
                                end
                            end
                        end
                    end
                end
            end
            for sex in 1:2
                for s_m in 1:NTypes_m
                    for s_f in 1:NTypes_f
                        for h_m in 1:hbin_m
                            for h_f in 1:hbin_f
                                for x_m in 1:xbin_m
                                    Threads.@threads for x_f in 1:xbin_f
                                        lock(lk) do
                                            FV_married[x_m,x_f,h_m,h_f,s_m,s_f,sex,t]=sum(p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t].*(vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,sex,t] .- log.(p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t])))
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            ###########################################
            # Expected values at beginning of stage 1
            ###########################################
            # Expected value for people who are MARRIED at the beginning of stage 1
            for sex in 1:2
                for s in 1:NTypes
                    for s_sp_orig in 1:NTypes_sp_orig
                        for h in 1:hbin
                            for h_sp_orig in 1:hbin_sp_orig
                                for phi in 1:phibin 
                                    x=phistate_list[phi,1]
                                    for phi_sp_orig in 1:phibin_sp_orig
                                        denom=0
                                        Threads.@threads for s_sp_dest in 1:NTypes
                                        # for s_sp_dest in 1:NTypes
                                            for h_sp_dest in 1:hbin
                                                # h_sp_dest_=hval[h_sp_dest]
                                                for phi_sp_dest in 1:phibin
                                                    x_sp_dest=phistate_list[phi_sp_dest,1]
                                                    if sex==1 
                                                        FV_married_= FV_married[x,x_sp_dest,h,h_sp_dest,s,s_sp_dest,sex,t]
                                                    else
                                                        FV_married_= FV_married[x_sp_dest,x,h_sp_dest,h,s_sp_dest,s,sex,t]
                                                    end
                                                    lock(lk) do
                                                        u_marry_marriedb4[phi,phi_sp_orig,phi_sp_dest,h,h_sp_orig,h_sp_dest,s,s_sp_orig,s_sp_dest,sex,t]=kappa(phi,phi_sp_orig,phi_sp_dest,h,h_sp_orig,h_sp_dest,s,s_sp_orig,s_sp_dest,sex) + FV_married_
                                                        denom=denom+exp(u_marry_marriedb4[phi,phi_sp_orig,phi_sp_dest,h,h_sp_orig,h_sp_dest,s,s_sp_orig,s_sp_dest,sex,t] )
                                                    end 
                                                end
                                            end
                                        end
                                        # Add the value of being single
                                        lock(lk) do
                                            u_single_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t]=kappa(phi,phi_sp_orig,0,h,h_sp_orig,0,s,s_sp_orig,0,sex) + FV[x,h, s, sex,t]
                                        end
                                        denom=denom+exp(u_single_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t])
                                        lock(lk) do
                                            FU_married[phi,phi_sp_orig, h, h_sp_orig, s, s_sp_orig, sex, t]= Beta* log(denom)
                                            p_st1_marriedb4[phi,phi_sp_orig,:,h,h_sp_orig,:,s,s_sp_orig,:,sex,t].=exp.(u_marry_marriedb4[phi,phi_sp_orig,:,h,h_sp_orig,:,s,s_sp_orig,:,sex,t])./denom
                                            p0_st1_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t]=exp(u_single_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t])/denom
                                        end
                                        # check=sum(p_st1_marriedb4[phi,phi_sp_orig,:,h,h_sp_orig,:,s,s_sp_orig,:,sex,t])+p0_st1_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t]
                                        # if abs(check-1)>1e-15
                                        #     # println("CCPs don't sum up to 1! $check")
                                        # end
                                    end
                                end
                            end
                        end
                    end
                end
            end

            # Expected value for people who are SINGLE at the beginning of stage 1
            for sex in 1:2
                for s in 1:NTypes
                    for h in 1:hbin
                        Threads.@threads for phi in 1:phibin 
                            x=phistate_list[phi,1]
                            denom=0
                            for s_sp_dest in 1:NTypes
                                for h_sp_dest in 1:hbin
                                    for phi_sp_dest in 1:phibin
                                        x_sp_dest=phistate_list[phi_sp_dest,1]
                                        if sex==1 
                                            FV_married_= FV_married[x,x_sp_dest,h,h_sp_dest,s,s_sp_dest,sex,t]
                                        else
                                            FV_married_= FV_married[x_sp_dest,x,h_sp_dest,h,s_sp_dest,s,sex,t]
                                        end
                                        lock(lk) do
                                            u_marry_single[phi,phi_sp_dest,h,h_sp_dest,s,s_sp_dest,sex,t]=kappa(phi,0,phi_sp_dest,h,0,h_sp_dest,s,0,s_sp_dest,sex) + FV_married_
                                        end
                                        denom=denom+exp( u_marry_single[phi,phi_sp_dest,h,h_sp_dest,s,s_sp_dest,sex,t])
                                    end
                                end
                            end
                            # Add the value of being single
                            lock(lk) do
                                u_single_single[phi,h,s,sex,t]=kappa(phi,0,0,h,0,0,s,0,0,sex) + FV[x,h, s, sex, t]
                            end
                            denom=denom + exp(u_single_single[phi,h,s,sex,t])
                            lock(lk) do
                                FU[phi, h, s, sex, t]= Beta* log(denom)
                                p_st1_singleb4[phi,:,h,:,s,:,sex,t].=exp.(u_marry_single[phi,:,h,:,s,:,sex,t])./denom
                                p0_st1_singleb4[phi,h,s,sex,t]=exp(u_single_single[phi,h,s,sex,t])/denom
                            end
                            check=sum(p_st1_singleb4[phi,:,h,:,s,:,sex,t])+p0_st1_singleb4[phi,h,s,sex,t]
                            if abs(check-1)>1e-15
                                # println("CCPs don't sum up to 1! $check")
                            end
                        end
                    end
                end
            end
            t -= 1
        end
        toc=time()
        duration=toc-tic

        CCP_true_sim=CCP_True_sim(
            p_st1_marriedb4,
            p0_st1_marriedb4,
            p_st1_singleb4,
            p0_st1_singleb4,
            p_st2,
            p_st2_cpl
        )
   
        p_st1_m_marriedb4=permutedims(p_st1_marriedb4[:,:,:,:,:,:,:,:,:,1,1:T2],(3,6,9,1,4,7,2,5,8,10))
        p_st1_f_marriedb4=permutedims(p_st1_marriedb4[:,:,:,:,:,:,:,:,:,2,1:T2],(3,6,9,1,4,7,2,5,8,10))
        p0_st1_m_marriedb4=permutedims(p0_st1_marriedb4[:,:,:,:,:,:,1,1:T2],(1,3,5,2,4,6,7))
        p0_st1_f_marriedb4=permutedims(p0_st1_marriedb4[:,:,:,:,:,:,2,1:T2],(1,3,5,2,4,6,7))
        p_st1_m_singleb4=permutedims(p_st1_singleb4[:,:,:,:,:,:,1,1:T2],(2,4,6,1,3,5,7))
        p_st1_f_singleb4=permutedims(p_st1_singleb4[:,:,:,:,:,:,2,1:T2],(2,4,6,1,3,5,7))
        p0_st1_m_singleb4=p0_st1_singleb4[:,:,:,1,1:T2]
        p0_st1_f_singleb4=p0_st1_singleb4[:,:,:,2,1:T2]
        p_st2_m=p_st2[:,:,:,:,1,1:T2]
        p_st2_f=p_st2[:,:,:,:,2,1:T2]
        p0_st2_m=p_st2_m[:,:,:,1,1:T2]
        p0_st2_f=p_st2_f[:,:,:,1,1:T2]
        p0_st2_cpl=permutedims(p_st2_cpl[:,:,:,:,:,:,1,1,1:T2],(1,3,5,2,4,6,7))
        p_st2_cpl=permutedims(p_st2_cpl[:,:,:,:,:,:,:,:,1:T2],(1,3,5,2,4,6,7,8,9))


        CCP_true=CCP_True(
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

        CCP_true=compress_structure(CCP_true, precision)
        CCP_true_sim=compress_structure(CCP_true_sim, precision)


        return CCP_true_sim, CCP_true
    end
    
    
    mutable struct CCP_True
        p0_st1_m_marriedb4:: Array
        p0_st1_f_marriedb4:: Array
        p0_st1_m_singleb4:: Array
        p0_st1_f_singleb4:: Array
        p0_st2_m:: Array
        p0_st2_f:: Array
        p0_st2_cpl:: Array
        p_st1_m_marriedb4:: Array
        p_st1_f_marriedb4:: Array
        p_st1_m_singleb4:: Array
        p_st1_f_singleb4:: Array
        p_st2_m:: Array
        p_st2_f:: Array
        p_st2_cpl:: Array
    end

    mutable struct CCP_True_sim
        p_st1_marriedb4:: Array
        p0_st1_marriedb4:: Array
        p_st1_singleb4:: Array
        p0_st1_singleb4:: Array
        p_st2:: Array
        p_st2_cpl:: Array
    end