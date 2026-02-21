    function f_data_genbus4_simulate(CCP_true_sim,N_)
        # In SOEP, people enter the panel when they cohabit with someone who is already in the panel. I need to replicate this in the data generation
        # Each individual starts as single and I generate their history. At each period, there is an exogenous probability that an individual is sampled by SOEP, as a household.
        # If in that period the individual is married, then their spouse is now considered an individual too and I generate a history for them too.
        # After being sampled, at each period, there is an exogenous probability that and individual will drop out of the panel, for whatever reason.
        ID_1 = [(n, 0) for n in 1:N_]
        # ID_1= collect(1:N_)
        Sex_1 = (rand(N_) .> P_male) .+1
        # State is the unobserved type
        State_1 = Int.(rand(N_) .> Pi_m).* (Sex_1.==1) .+ Int.(rand(N_) .> Pi_f).* (Sex_1.==2) .+1
        State_sp_orig_1=Array{Union{Missing, Int64}}(undef, N_, T+1)
        State_sp_dest_1=Array{Union{Missing, Int64}}(undef, N_, T)
        Y_1 = Array{Union{Missing, Int64}}(undef, N_, T)
        Y_sp_dest_1=Array{Union{Missing, Int64}}(undef, N_, T)
        Y_sp_orig_1=Array{Union{Missing, Int64}}(undef, N_, T)
        # Hstate is the randomly drawn position index on the grid of values of h
        Hstate_1 = Int.(ceil.(length(hval) * rand(N_)))
        Hstate_sp_orig_1=Array{Union{Missing, Int64}}(undef, N_, T+1)
        Hstate_sp_dest_1=Array{Union{Missing, Int64}}(undef, N_, T+1)
        H_1 = hval[Hstate_1]
        H_sp_dest_1=Array{Union{Missing, precision}}(undef, N_, T+1)
        H_sp_orig_1=Array{Union{Missing, precision}}(undef, N_, T+1)
        Phi_1 = Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phi_1[:,1,1].=0
        Phi_sp_orig_1=Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phi_sp_dest_1=Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phistate_1 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Phistate_1[:,1].=1
        Phistate_sp_orig_1 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Phistate_sp_dest_1 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Married_1=Array{Union{Missing, Int64}}(undef, N_, T)
        Sampled_1=zeros(Int,N_,T )
        Draw_1_stage2 = rand(N_, T )
        Draw_1_transition = rand(N_, T )
        Draw_1_stage1 = rand(N_, T )
        Draw_1_sample = rand(N_, T )
        Draw_1_drop = rand(N_, T )
        Temp_1=zeros(N_,T)

        ID_2=Array{Union{Missing, Tuple}}(missing, N_)
        Sex_2 = Array{Union{Missing, Int64}}(undef, N_)
        # State is the unobserved type
        State_2 = Array{Union{Missing, Int64}}(undef, N_)
        State_sp_orig_2=Array{Union{Missing, Int64}}(undef, N_, T+1)
        State_sp_dest_2=Array{Union{Missing, Int64}}(undef, N_, T)
        Y_2 = Array{Union{Missing, Int64}}(undef, N_, T)
        Y_sp_dest_2= Array{Union{Missing, Int64}}(undef, N_, T)
        Y_sp_orig_2= Array{Union{Missing, Int64}}(undef, N_, T)
        # Hstate is the randomly drawn position index on the grid of values of h
        Hstate_2 = Array{Union{Missing, Int64}}(undef, N_)
        Hstate_sp_orig_2=Array{Union{Missing, Int64}}(undef, N_, T+1)
        Hstate_sp_dest_2=Array{Union{Missing, Int64}}(undef, N_, T+1)
        H_2 = Array{Union{Missing, precision}}(undef, N_)
        H_sp_dest_2=Array{Union{Missing, precision}}(undef, N_, T+1)
        H_sp_orig_2=Array{Union{Missing, precision}}(undef, N_, T+1)
        Phi_2 = Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phi_sp_orig_2=Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phi_sp_dest_2=Array{Union{Missing, precision}}(undef, N_, T+1, phidim)
        Phistate_2 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Phistate_sp_orig_2 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Phistate_sp_dest_2 = Array{Union{Missing, Int64}}(undef, N_, T+1)
        Married_2=Array{Union{Missing, Int64}}(undef, N_, T)
        Sampled_2=zeros(Int,N_,T )
        Draw_2_stage2 = rand(precision,(N_, T) )
        Draw_2_transition = rand(precision,(N_, T) )
        Draw_2_stage1 = rand(precision,(N_, T) )
        Draw_2_sample = rand(precision, (N_, T) )
        Draw_2_drop = rand(precision, (N_, T) )
        Temp_1=zeros(precision, (N_,T))

        temp1=reshape(collect(1:phibin*hbin*NTypes),(phibin,hbin,NTypes))
        temp2=reshape(collect(1:(J+1)^2),(J+1,J+1))
        temp3=reshape(collect(1:phibin_m*phibin_f),(phibin_m,phibin_f))

        for n in 1:N_
            # println(n)
            # Convert position index on H space into H value
            h=Hstate_1[n]
            h_ = hval[h]
            H_1[n] = h_
            sex=Sex_1[n]
            for t in 1:T 
                # println(t)
                s=State_1[n]
                phi=Phistate_1[n,t]
                x=phistate_list[phi,1]
                if t==1
                    married_tminus1=0
                else
                    married_tminus1=Married_1[n,t-1]
                    if married_tminus1==1
                        if COND_ORIG_SPOUSE
                            s_sp_orig=Int(State_sp_orig_1[n,t])
                            h_sp_orig=Hstate_sp_orig_1[n,t]
                            phi_sp_orig=Phistate_sp_orig_1[n,t]
                        else
                            s_sp_orig=1
                            h_sp_orig=1
                            phi_sp_orig=1
                        end
                    end
                end

                ##########################
                # Stage 1: Choose partner
                ##########################
                if married_tminus1==1
                    psi_marry=CCP_true_sim.p_st1_marriedb4[phi,phi_sp_orig,:,h,h_sp_orig,:,s,s_sp_orig,:,sex,t]
                    psi_single=CCP_true_sim.p0_st1_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t]
                else
                    psi_marry=CCP_true_sim.p_st1_singleb4[phi,:,h,:,s,:,sex,t]
                    psi_single=CCP_true_sim.p0_st1_singleb4[phi,h,s,sex,t]
                end
                cumulprob=cumsum(vcat(psi_single,psi_marry[:]))
                spouse=searchsortedfirst(cumulprob,Draw_1_stage1[n, t])
                spouse=min(spouse,length(cumulprob))
                married=spouse>1

                Married_1[n,t]=married
                if married==1
                    # Because of how I did vcat above, spouse=1 means the person remains single
                    spouse=spouse-1
                    # Spouse now indicates the new spouse position in the array of size (phibin,hbin,NTypes)
                    Phistate_sp_dest_1[n,t]=findall(spouse.==temp1)[1][1]
                    Phi_sp_dest_1[n,t,:]=phival[Phistate_sp_dest_1[n,t],:]
                    Hstate_sp_dest_1[n,t]=findall(spouse.==temp1)[1][2]
                    H_sp_dest_1[n,t]=hval[Hstate_sp_dest_1[n,t]]
                    State_sp_dest_1[n,t]=findall(spouse.==temp1)[1][3]
                end

                ##########################
                # Stage 1.1: sampling
                ##########################
                if t==1
                    Sampled_1[n,t]=(rand()<P_sample)*1
                elseif Sampled_1[n,t-1]==0 
                    Sampled_1[n,t]=(rand()<P_sample)*1
                elseif Sampled_1[n,t-1]==1
                    Sampled_1[n,t]=1
                end

                ##########################
                # Stage 2: Other choices
                ##########################
                if married==1
                    if sex==1
                        phi_m=phi
                        phi_f=Phistate_sp_dest_1[n,t]
                        h_m=h
                        h_f=Hstate_sp_dest_1[n,t]
                        s_m=s
                        s_f=State_sp_dest_1[n,t]
                    else
                        phi_f=phi
                        phi_m=Phistate_sp_dest_1[n,t]
                        h_f=h
                        h_m=Hstate_sp_dest_1[n,t]
                        s_f=s
                        s_m=State_sp_dest_1[n,t]
                    end
                    x_m=phistate_list[phi_m,1]
                    x_f=phistate_list[phi_f,1]
                    p=CCP_true_sim.p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t] #Remember that in this array, j=1 means replacing the engine, it's the normalized action
                    temp=searchsortedfirst(cumsum(p[:]),Draw_1_stage2[n, t])
                    temp=min(temp,length(cumsum(p[:]))) #When operating at low precisions, CCPs sometimes add up to slightly below 1 and the code fails if a draw is above that sum
                    j_m=findall(temp.==temp2)[1][1]
                    j_f=findall(temp.==temp2)[1][2]
                    if sex==1
                        Y_1[n,t]=j_m
                        Y_sp_dest_1[n,t]=j_f
                        if t>1
                            Y_sp_orig_1[n,t-1]=j_f
                        end
                    else
                        Y_1[n,t]=j_f
                        Y_sp_dest_1[n,t]=j_m
                        if t>1
                            Y_sp_orig_1[n,t-1]=j_m
                        end
                    end
                else  
                    p=CCP_true_sim.p_st2[x,h,s,:,sex,t]
                    j=searchsortedfirst(cumsum(p[:]),Draw_1_stage2[n, t])
                    j=min(j,length(cumsum(p[:])))
                    Y_1[n,t]=j
                end

                ###################################
                # State transition to next period
                ###################################
                if married==1
                    temp=searchsortedfirst(xtranc_married[:,:,x_m,x_f,h_m,h_f,j_m,j_f][:],Draw_1_transition[n,t])
                    temp=min(temp,length(xtranc_married[:,:,x_m,x_f,h_m,h_f,j_m,j_f][:])) #When operating at low precisions, CCPs sometimes add up to slightly below 1 and the code fails if a draw is above that sum 
                    Temp_1[n,t]=temp
                    position = findfirst(x -> x == temp, temp3)
                    phi_m=position[1]
                    phi_f=position[2]
                    if sex==1
                        Phistate_1[n, t+1] = phi_m
                        Phi_1[n,t+1,:]=phival[phi_m,:]
                        Phistate_sp_orig_1[n, t+1] = phi_f
                        Phi_sp_orig_1[n,t+1,:]=phival[phi_f,:]
                    else
                        Phistate_1[n, t+1] = phi_f
                        Phi_1[n,t+1,:]=phival[phi_f,:]
                        Phistate_sp_orig_1[n, t+1] = phi_m
                        Phi_sp_orig_1[n,t+1,:]=phival[phi_m,:]
                    end
                        State_sp_orig_1[n,t+1]=State_sp_dest_1[n,t]
                        Hstate_sp_orig_1[n,t+1]=Hstate_sp_dest_1[n,t]
                        H_sp_orig_1[n,t+1]=H_sp_dest_1[n,t]
                else
                    temp=searchsortedfirst(xtranc[x,:,h,j],Draw_1_transition[n,t])
                    temp=min(temp,length(xtranc[x,:,h,j]))
                    Phistate_1[n, t+1]=temp
                    Phi_1[n, t+1,:] = phival[Phistate_1[n, t+1],:]
                end
                Dropped=rand()<P_drop
                if Dropped==1
                    break
                end
            end

            #############################################################################################################################
            # Now, I simulate the history of their spouses at the time they were first sampled
            #############################################################################################################################
            
            t_sample=findfirst(x -> x == 1, Sampled_1[n,:])
            if isnothing(t_sample)
                continue
            end
            if Married_1[n,t_sample]==1
                ID_2[n]=(-ID_1[n][1],0)
                # Convert position index on H space into H value
                h=Hstate_sp_dest_1[n,t_sample]
                h_ = hval[h]
                H_2[n] = h_
                Hstate_2[n]=h
                Sex_2[n]=(Sex_1[n]==1)*2+(Sex_1[n]==2)
                sex=Int(Sex_2[n])
                s=State_sp_dest_1[n,t_sample]
                State_2[n]=s
                phi_=Phi_sp_dest_1[n, t_sample,:]
                Phi_2[n,t_sample,:]=phi_
                phi=Phistate_sp_dest_1[n,t_sample]
                Phistate_2[n,t_sample]=phi
                Married_2[n,t_sample]=1
                married=1
                Phistate_sp_dest_2[n,t_sample]=Phistate_1[n,t_sample]
                Phi_sp_dest_2[n,t_sample,:]=phival[ Phistate_sp_dest_2[n,t_sample],:]
                Hstate_sp_dest_2[n,t_sample]=Hstate_1[n]
                H_sp_dest_2[n,t_sample]=hval[Hstate_sp_dest_2[n,t_sample]]
                State_sp_dest_2[n,t_sample]= State_1[n]
                for t in t_sample:T 
                    Sampled_2[n,t]=1
                    if t>t_sample
                        phi=Phistate_2[n,t]
                        x=phistate_list[phi,1]
                        married_tminus1=Married_2[n,t-1]
                        if married_tminus1==1
                            if COND_ORIG_SPOUSE
                                s_sp_orig=Int(State_sp_orig_2[n,t])
                                h_sp_orig=Hstate_sp_orig_2[n,t]
                                phi_sp_orig=Phistate_sp_orig_2[n,t]
                            else
                                s_sp_orig=1
                                h_sp_orig=1
                                phi_sp_orig=1
                            end
                        end
                    
                        ##########################
                        # Stage 1: Choose partner
                        ##########################
                        if married_tminus1==1
                            psi_marry=CCP_true_sim.p_st1_marriedb4[phi,phi_sp_orig,:,h,h_sp_orig,:,s,s_sp_orig,:,sex,t]
                            psi_single=CCP_true_sim.p0_st1_marriedb4[phi,phi_sp_orig,h,h_sp_orig,s,s_sp_orig,sex,t]
                        else
                            psi_marry=CCP_true_sim.p_st1_singleb4[phi,:,h,:,s,:,sex,t]
                            psi_single=CCP_true_sim.p0_st1_singleb4[phi,h,s,sex,t]
                        end
                        cumulprob=cumsum(vcat(psi_single,psi_marry[:]))
                        spouse=searchsortedfirst(cumulprob,Draw_2_stage1[n, t])
                        spouse=min(spouse,length(cumulprob))
                        married=spouse>1

                        Married_2[n,t]=married
                        if married==1
                            # Because of how I did vcat above, spouse=1 means the person remains single
                            spouse=spouse-1
                            # Spouse now indicates the new spouse position in the array of size (xbin,hbin,NTypes)
                            Phistate_sp_dest_2[n,t]=findall(spouse.==temp1)[1][1]
                            Phi_sp_dest_2[n,t,:]=phival[Phistate_sp_dest_2[n,t],:]
                            Hstate_sp_dest_2[n,t]=findall(spouse.==temp1)[1][2]
                            H_sp_dest_2[n,t]=hval[Hstate_sp_dest_2[n,t]]
                            State_sp_dest_2[n,t]=findall(spouse.==temp1)[1][3]
                        end

                    end

                    ##########################
                    # Stage 2: Other choices
                    ##########################
                    if married==1
                        if sex==1
                            phi_m=phi
                            phi_f=Phistate_sp_dest_2[n,t]
                            h_m=h
                            h_f=Hstate_sp_dest_2[n,t]
                            s_m=s
                            s_f=State_sp_dest_2[n,t]
                        else
                            phi_f=phi
                            phi_m=Phistate_sp_dest_2[n,t]
                            h_f=h
                            h_m=Hstate_sp_dest_2[n,t]
                            s_f=s
                            s_m=State_sp_dest_2[n,t]
                        end
                        x_m=phistate_list[phi_m,1]
                        x_f=phistate_list[phi_f,1]
                        p=CCP_true_sim.p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,:,:,t] #Remember that in this array, j=1 means replacing the engine, it's the normalized action
                        temp=searchsortedfirst(cumsum(p[:]),Draw_2_stage2[n, t])
                        temp=min(temp,length(cumsum(p[:])))
                        j_m=findall(temp.==temp2)[1][1]
                        j_f=findall(temp.==temp2)[1][2]
                        if sex==1
                            Y_2[n,t]=j_m
                            Y_sp_dest_2[n,t]=j_f
                            if t>1
                                Y_sp_orig_2[n,t-1]=j_f
                            end
                        else
                            Y_2[n,t]=j_f
                            Y_sp_dest_2[n,t]=j_m
                            if t>1
                                Y_sp_orig_2[n,t-1]=j_m
                            end
                        end
                    else  
                        p=CCP_true_sim.p_st2[x,h,s,:,sex,t]
                        j=searchsortedfirst(cumsum(p[:]),Draw_2_stage2[n, t])
                        j=min(j,length(cumsum(p[:])))
                        Y_2[n,t]=j
                    end
                    

                    ###################################
                    # State transition to next period
                    ###################################
                    if married==1
                        # j_m1=j_m*(j_m>0)+(j_m==0)
                        # phi_m1=phi_m*(j_m>0)+(j_m==0)
                        # j_f1=j_f*(j_f>0)+(j_f==0)
                        # phi_f1=phi_f*(j_f>0)+(j_f==0)
                        # temp=searchsortedfirst(xtranc_married[:,:,phi_m1,phi_f1,h_m,h_f,j_m1,j_f1][:],Draw_2_transition[n,t])
                        temp=searchsortedfirst(xtranc_married[:,:,x_m,x_f,h_m,h_f,j_m,j_f][:],Draw_2_transition[n,t])
                        temp=min(temp,length(xtranc_married[:,:,x_m,x_f,h_m,h_f,j_m,j_f][:])) #When operating at low precisions, CCPs sometimes add up to slightly below 1 and the code fails if a draw is above that sum
                        # phi_m=mod(temp,xbin)*(mod(temp,xbin)>0)+xbin*(mod(temp,xbin)==0)
                        # phi_f=Int(ceil(temp/xbin))
                        position = findfirst(x -> x == temp, temp3)
                        phi_m=position[1]
                        phi_f=position[2]
                        if sex==1
                            Phistate_2[n, t+1] = phi_m
                            Phi_2[n,t+1,:]=phival[phi_m,:]
                            Phistate_sp_orig_2[n, t+1] = phi_f
                            Phi_sp_orig_2[n,t+1,:]=phistate_list[phi_f,:]
                        else
                            Phistate_2[n, t+1] = phi_f
                            Phi_2[n,t+1,:]=phival[phi_f,:]
                            Phistate_sp_orig_2[n, t+1] = phi_m
                            Phi_sp_orig_2[n,t+1,:]=phival[phi_m,:]
                        end
                        State_sp_orig_2[n,t+1]=State_sp_dest_2[n,t]
                        Hstate_sp_orig_2[n,t+1]=Hstate_sp_dest_2[n,t]
                        H_sp_orig_2[n,t+1]=H_sp_dest_2[n,t]
                    else
                        temp=searchsortedfirst(xtranc[x,:,h,j],Draw_2_transition[n,t])
                        temp=min(temp,length(xtranc[x,:,h,j]))
                        Phistate_2[n, t+1]=temp
                        Phi_2[n, t+1,:] = phival[Phistate_2[n, t+1],:]
                    end

                    Dropped=rand()<P_drop
                    if Dropped==1
                        break
                    end

                end
            end
        end

        # Keep only individuals that were sampled at all
        # sampled_1 = findall(row -> any(x -> x == 1, row), eachrow(Sampled_1))

        #keep only individuals for which we have at least T3 periods
        sampled_1 = findall(row -> sum(row)>=T3, eachrow(Sampled_1))
        sampled_2 = maximum(Sampled_2[sampled_1,:], dims=2)[:]

        Sampled_1=Sampled_1[sampled_1,1:T2]
        Sampled_2=Sampled_2[sampled_1,1:T2]
     
        State_1=State_1[sampled_1]
        State_2=State_2[sampled_1]
        State_2[sampled_2.==0].=missing
        State_sp_dest_1=State_sp_dest_1[sampled_1,1:T2]
        State_sp_dest_2=State_sp_dest_2[sampled_1,1:T2]
        State_sp_dest_2[sampled_2.==0,:].=missing
        FirstSampled_1= [findfirst(x -> x == 1, row) for row in eachrow(Sampled_1)][:]
        Y_1=Y_1[sampled_1,1:T2]
        Y_1[Sampled_1.==0].=missing
        Y_sp_dest_1=Y_sp_dest_1[sampled_1,1:T2]
        Y_sp_dest_1[Sampled_1.==0].=missing
        Y_2=Y_2[sampled_1,1:T2]
        Y_2[Sampled_2.==0].=missing
        Y_sp_dest_2=Y_sp_dest_2[sampled_1,1:T2]
        Y_sp_dest_2[Sampled_2.==0].=missing
        Phi_1=Phi_1[sampled_1,1:T2,:]
        Phi_1[repeat(Sampled_1,1,1,phidim).==0].=missing
        Phi_sp_dest_1=Phi_sp_dest_1[sampled_1,1:T2,:]
        Phi_sp_dest_1[repeat(Sampled_1,1,1,phidim).==0].=missing
        Phi_sp_orig_1=Phi_sp_orig_1[sampled_1,1:T2,:]
        Phi_sp_orig_1[repeat(Sampled_1,1,1,phidim).==0].=missing
        Phi_2=Phi_2[sampled_1,1:T2,:]
        Phi_2[repeat(Sampled_2,1,1,phidim).==0].=missing
        Phi_sp_dest_2=Phi_sp_dest_2[sampled_1,1:T2,:]
        Phi_sp_dest_2[repeat(Sampled_2,1,1,phidim).==0].=missing
        Phi_sp_orig_2=Phi_sp_orig_2[sampled_1,1:T2,:]
        Phi_sp_orig_2[repeat(Sampled_2,1,1,phidim).==0].=missing
        Phistate_1=Phistate_1[sampled_1,1:T2]
        Phistate_1[Sampled_1.==0].=missing
        Phistate_sp_dest_1=Phistate_sp_dest_1[sampled_1,1:T2]
        Phistate_sp_dest_1[Sampled_1.==0].=missing
        Phistate_sp_orig_1=Phistate_sp_orig_1[sampled_1,1:T2]
        Phistate_sp_orig_1[Sampled_1.==0].=missing
        Phistate_2=Phistate_2[sampled_1,1:T2]
        Phistate_2[Sampled_2.==0].=missing
        Phistate_sp_dest_2=Phistate_sp_dest_2[sampled_1,1:T2]
        Phistate_sp_dest_2[Sampled_2.==0].=missing
        Phistate_sp_orig_2=Phistate_sp_orig_2[sampled_1,1:T2]
        Phistate_sp_orig_2[Sampled_2.==0].=missing
        H_1=H_1[sampled_1]
        H_sp_dest_1=H_sp_dest_1[sampled_1,1:T2]
        H_sp_dest_1[Sampled_1.==0].=missing
        H_sp_orig_1=H_sp_orig_1[sampled_1,1:T2]
        H_sp_orig_1[Sampled_1.==0].=missing
        H_2=H_2[sampled_1]
        H_2[sampled_2.==0].=missing
        H_sp_dest_2=H_sp_dest_2[sampled_1,1:T2]
        H_sp_dest_2[Sampled_2.==0].=missing
        H_sp_orig_2=H_sp_orig_2[sampled_1,1:T2]
        H_sp_orig_2[Sampled_2.==0].=missing
        Hstate_1=Hstate_1[sampled_1]
        Hstate_sp_dest_1=Hstate_sp_dest_1[sampled_1,1:T2]
        Hstate_sp_dest_1[Sampled_1.==0].=missing
        Hstate_sp_orig_1=Hstate_sp_orig_1[sampled_1,1:T2]
        Hstate_sp_orig_1[Sampled_1.==0].=missing
        Hstate_2=Hstate_2[sampled_1]
        Hstate_2[sampled_2.==0].=missing
        Hstate_sp_dest_2=Hstate_sp_dest_2[sampled_1,1:T2]
        Hstate_sp_dest_2[Sampled_2.==0].=missing
        Hstate_sp_orig_2=Hstate_sp_orig_2[sampled_1,1:T2]
        Hstate_sp_orig_2[Sampled_2.==0].=missing
        Sex_1=Sex_1[sampled_1]
        Sex_2=Sex_2[sampled_1]
        Sex_2[sampled_2.==0].=missing
        ID_1=ID_1[sampled_1]
        ID_2=ID_2[sampled_1]
        ID_2[sampled_2.==0].=missing
        Married_1=Married_1[sampled_1,1:T2]
        Married_1[Sampled_1.==0].=missing
        Married_2=Married_2[sampled_1,1:T2]
        Married_2[Sampled_2.==0].=missing


        N_=size(Sampled_1)[1]


        # Additionally, all _sp_orig variables should have missing values in the first period an agent is in the sample
        for n in 1:N_
            tfirst=FirstSampled_1[n]
            Phi_sp_orig_1[n,tfirst,:].=missing
            Phi_sp_orig_2[n,tfirst,:].=missing
            Phistate_sp_orig_1[n,tfirst]=missing
            Phistate_sp_orig_2[n,tfirst]=missing
            H_sp_orig_1[n,tfirst]=missing
            H_sp_orig_2[n,tfirst]=missing
            Hstate_sp_orig_1[n,tfirst]=missing
            Hstate_sp_orig_2[n,tfirst]=missing
        end
    
    
        # Create ID for spouses too and mapping from individuals to spouse for each t
        ID_sp_dest_1=Array{Union{Missing, Tuple}}(missing, N_, T2)
        ID_sp_dest_2=Array{Union{Missing, Tuple}}(missing, N_, T2)
        
        for n in 1:N_
            # println(n)
            id_sp_1=ifelse(ismissing(ID_2[n]),ID_1[n],ID_2[n])
            id_sp_2=ID_1[n]
            condition=false
            for t in 1:T2
                # println(t)
                if t<FirstSampled_1[n]
                    continue
                end
                if (t==FirstSampled_1[n]) & (Married_1[n,t]===1)
                    ID_sp_dest_1[n,t]=id_sp_1
                    ID_sp_dest_2[n,t]=id_sp_2
                else
                    # condition_1 and condition_2 indicate whether agent 1 or 2 respectively changed partner since the previous period
                    condition_1=(Phistate_sp_orig_1[n,t]!==Phistate_sp_dest_1[n,t]) | (Hstate_sp_orig_1[n,t]!==Hstate_sp_dest_1[n,t]) | (State_sp_orig_1[n,t]!==State_sp_dest_1[n,t])
                    condition_2=(Phistate_sp_orig_2[n,t]!==Phistate_sp_dest_2[n,t]) | (Hstate_sp_orig_2[n,t]!==Hstate_sp_dest_2[n,t]) | (State_sp_orig_2[n,t]!==State_sp_dest_2[n,t])
                    
                    # condition indicates that the couple formed of agent 1 and agent 2 at sampling broke up, because either of the two agents matched to somebody else.
                    if t>1
                        condition=(ID_sp_dest_1[n,t-1]===ID_2[n]) & (ID_sp_dest_1[n,t-1]!==missing) & (condition_1 | condition_2)
                    end
                    # If the original couple broke up, initialize the spouse ID as (i,0) or (-i,0)
                    if condition==true
                        id_sp_1=(ID_1[n][1],0)
                        id_sp_2=(ID_2[n][1],0)
                    end
                    # If agent 1 is married at t and either changed partner since t-1 or agent 2 broke up with 1 since t-1, then increase the spouse counter by 1
                    if (condition_1 | condition) & (Married_1[n,t]===1)
                        id_sp_1=(id_sp_1[1],id_sp_1[2]+1)
                    end
                    # same thing for agent 2
                    if (condition_2 | condition) & (Married_2[n,t]===1)
                        id_sp_2=(id_sp_2[1],id_sp_2[2]+1)
                    end
                    
                    if (Married_1[n,t]===1) & (Sampled_1[n,t].==1)
                        #  If agent 1 is in sample and married, paste id_sp_1, however defined at this point, into ID_sp_dest_1
                        ID_sp_dest_1[n,t]=id_sp_1
                    elseif (Married_1[n,t]===0) & (Sampled_1[n,t].==1)
                        # If agent 1 is in sample and single, paste their own id (ID_1[n][1],0), into ID_sp_dest_1
                        ID_sp_dest_1[n,t]=(ID_1[n][1],0)
                    elseif (Sampled_1[n,t].==0)
                        # if agent 1 is not in sample at this t, ID_sp_dest_1 should be missing
                        ID_sp_dest_1[n,t]=missing
                    end
                    # same thing for agent 2
                    if (Married_2[n,t]===1) & (Sampled_2[n,t].==1)
                        ID_sp_dest_2[n,t]=id_sp_2
                    elseif (Married_2[n,t]===0) & (Sampled_2[n,t].==1)
                        ID_sp_dest_2[n,t]=(ID_2[n][1],0)
                    elseif (Sampled_2[n,t].==0)
                        ID_sp_dest_2[n,t]=missing
                    end
                end
                # condition can be equal to true only once per history, so switch it back to false
                condition=false
            end
        end
    
        dt_mat=Data_Matrices(
        Y_1, 
        Y_sp_dest_1, 
        Y_sp_orig_1,
        Phi_1, 
        Phi_sp_dest_1, 
        Phi_sp_orig_1, 
        H_1, 
        H_sp_dest_1,
        H_sp_orig_1,
        Phistate_1, 
        Phistate_sp_dest_1, 
        Phistate_sp_orig_1, 
        Hstate_1, 
        Hstate_sp_dest_1, 
        Hstate_sp_orig_1, 
        State_1, 
        State_sp_dest_1, 
        State_sp_orig_1,
        ID_1, 
        Sampled_1, 
        Married_1, 
        Sex_1,
        Y_2, 
        Y_sp_dest_2, 
        Y_sp_orig_2,
        Phi_2, 
        Phi_sp_dest_2, 
        Phi_sp_orig_2, 
        H_2, 
        H_sp_dest_2, 
        H_sp_orig_2,
        Phistate_2, 
        Phistate_sp_dest_2, 
        Phistate_sp_orig_2, 
        Hstate_2, 
        Hstate_sp_dest_2, 
        Hstate_sp_orig_2, 
        State_2, 
        State_sp_dest_2, 
        State_sp_orig_2, 
        ID_2, 
        Sampled_2, 
        Married_2, 
        Sex_2, 
        Adj,
        ID_sp_dest_1,
        ID_sp_dest_2,
        FirstSampled_1)
        return dt_mat
    end

    
    struct Data_Matrices
        Y_1:: Array 
        Y_sp_dest_1:: Array 
        Y_sp_orig_1:: Array 
        Phi_1:: Array 
        Phi_sp_dest_1:: Array 
        Phi_sp_orig_1:: Array 
        H_1:: Array 
        H_sp_dest_1:: Array
        H_sp_orig_1:: Array
        Phistate_1:: Array 
        Phistate_sp_dest_1:: Array 
        Phistate_sp_orig_1:: Array 
        Hstate_1:: Array 
        Hstate_sp_dest_1:: Array 
        Hstate_sp_orig_1:: Array 
        State_1:: Array 
        State_sp_dest_1:: Array 
        State_sp_orig_1:: Array
        ID_1:: Array 
        Sampled_1:: Array 
        Married_1:: Array 
        Sex_1:: Array
        Y_2:: Array 
        Y_sp_dest_2:: Array 
        Y_sp_orig_2:: Array 
        Phi_2:: Array 
        Phi_sp_dest_2:: Array 
        Phi_sp_orig_2:: Array 
        H_2:: Array 
        H_sp_dest_2:: Array 
        H_sp_orig_2:: Array
        Phistate_2:: Array 
        Phistate_sp_dest_2:: Array 
        Phistate_sp_orig_2:: Array 
        Hstate_2:: Array 
        Hstate_sp_dest_2:: Array 
        Hstate_sp_orig_2:: Array 
        State_2:: Array 
        State_sp_dest_2:: Array 
        State_sp_orig_2:: Array 
        ID_2:: Array 
        Sampled_2:: Array 
        Married_2:: Array 
        Sex_2:: Array 
        Adj:: Array
        ID_sp_dest_1:: Array
        ID_sp_dest_2:: Array
        FirstSampled_1:: Array
    end