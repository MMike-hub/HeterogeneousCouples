# Bus data generation
# Adds another variable
# Change from orig folder to go with Jon's parameters for Z and X from 3_z

function genbus4(alpha, N, T, xtran, xtranc, xtran_married, xtranc_married, xbin, zbin, xval, zval, T2)
    """
    Generates decision data for a CCP Model
        alpha [vector float]: Vector of coefficients to parameterize the model 
        N [int]: Number of individuals 
        T [int]: Number of periods to simulate 
        xtran [Matrix float]: Transition probability matrix
        xval [int]: number of observed states
    
    Returns: 
        Y [matrix int]: NxT matrix of decisions of each individual at time t  
        X [matrix float]: NxT matrix of observed states of each individual at time t
        Xstate [matrix int]: NxT matrix of index of states of each individual at time t
        State [matrix int]: NxT matrix of unobs. state of each individual at time t
        FV [matrix float]: FOLLOWING OUR PAPER'S NOTATION. FV is the expected value at stage 2,  i.e. after the marriage decision is taken and before the other choice shocks are realized.
        FU [array float]: FOLLOWING OUR PAPER'S NOTATION. FU is the expected value at stage 1, i.e. before the marriage decision
    """
    global Adj, FVT

    # alpha=parameters
    # N=sample size
    # T=time horizon
    # xtran=transition matrix
    # xtranc=cumulative transition probability, also reshaped with respect to xtran
    # xbin, zbin=size of state space in each dimension
    # xval, zval= discretized state space dimensions
    Beta = alpha[5]
    # Pi=proportion of unobserved types s=1
    Pi = alpha[6]
    eul = 0.577215665
    # Adj is a set of year fixed effects
    Adj = zeros(T + 11)

    t = 2
    while t < T + 12
        Adj[t] = 0.7 * Adj[t-1] + 0.5 * randn()
        t += 1
    end

    # Backwards induction to solve model, i.e., find the continuation values big V
    # Notice that the lifespan is assumed to be T (column T+1 in FV is all zeros, i.e., utility after last period is zero)
    # FV= forward values. First dimension is destination state [(x_t+1,z_t+1)], second dimension is choice at time t, third dimension is time t+1
    FV = zeros(xbin,zbin, NTypes, NSex, T  + 1)
    FV_married=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes, NSex, T  + 1)
    # Due to transition costs in and out of marriage, the expected value at the beginning of stage 1 also depends on the state of the partner.
    # FU contains the expected value of single agents at the beginning of stage 1 
    FU = zeros(xbin, zbin, NTypes, NSex, T  + 1)
    # FU_married contains the expected value of married agents at the beginning of stage 1, as a function of who they are married to
    FU_married=zeros(xbin,xbin,zbin, zbin, NTypes, NTypes, NSex, T  + 1)
    vJ=zeros(xbin,zbin,NTypes,J,NSex,T)
    v0=zeros(xbin,zbin,NTypes,NSex,T)
    pJ_single=zeros(xbin,zbin,NTypes,J,NSex,T)
    p0_single=zeros(xbin,zbin,NTypes,NSex,T)
    v0J_married=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,J+1,J+1,NSex,T)
    # p0J_married gives the probability that a couple of type (x,x_sp,z,z_sp,s,s_sp) at time t takes a discrete decision (j,j_sp)
    p0J_married=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,J+1,J+1,T)
    # lambda contains the pareto-weight (bargaining power) of the male in each couple
    # lambda=ones(xbin,xbin,zbin,zbin,NTypes,NTypes,J+1).*0.5
    lambda=0.5
     # u_ is the conditional valuation function of the marriage choices at stage 1
     u_marry_married=zeros(xbin,xbin,xbin,zbin,zbin,zbin,NTypes,NTypes,NTypes,NSex,T) # The dimensions represent (x,x_spouse_origin,x_spouse_destination,z,z_spouse_origin,z_spouse_destination,s,s_spouse_origin,s_spouse_destination,sex)
     p_marry_married=zeros(xbin,xbin,xbin,zbin,zbin,zbin,NTypes,NTypes,NTypes,NSex,T)
     u_single_married=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,NSex,T) # The dimensions represent (x,x_spouse_origin,z,z_spouse_origin,s,s_spouse_origin,sex)
     p_single_married=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,NSex,T)
     u_marry_single=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,NSex,T) # The dimensions represent (x,x_spouse_destination,z,z_spouse_destination,s,s_spouse_destination,sex)
     p_marry_single=zeros(xbin,xbin,zbin,zbin,NTypes,NTypes,NSex,T)
     u_single_single=zeros(xbin,zbin,NTypes,NSex,T) # The dimensions represent (x,z,s,sex)
     p_single_single=zeros(xbin,zbin,NTypes,NSex,T) 
     
    # The kappa function computes transition costs.
    function kappa(x_,x_sp_orig_,x_sp_dest_,z_,z_sp_orig_,z_sp_dest_,s,s_sp_orig,s_sp_dest,sex)
        local kappa
        if (x_sp_dest_!=x_sp_orig_) || (z_sp_orig_!=z_sp_dest_) || (s_sp_orig!=s_sp_dest)
            kappa=1
        else
            kappa=0
        end
        return kappa
    end

    function util_fn(x,z,s,j,sex,t)
        local util_fn, x_, z_
        x_=xval[x]
        z_=zval[z]
        util_fn= alpha[1] + Adj[t] + alpha[2] * x_ + alpha[3] * (s==2) 
        return util_fn
    end

    function util_married_fn(x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,sex,t)
        local util_fn, x_m_, z_m_, x_f_, z_f_
        x_m_=xval[x_m]
        z_m_=zval[z_m]
        x_f_=xval[x_f]
        z_f_=zval[z_f]
        if sex==1 
            x_=x_m_
            s=s_m
            j=j_m
        else
            x_=x_f_
            s=s_f
            j=j_f
        end
        util_fn= alpha[1] + Adj[t] + alpha[2] * x_ + alpha[3] * (s==2) + alpha[4] * j
        return util_fn
    end


    tic=time()
    lk = ReentrantLock()
    t = T 
    while t >= 1
        #############################
        # Choices as single
        #############################
        for s in 1:2
            for z in 1:zbin
                Threads.@threads for x in 1:xbin 
                    for sex=1:2
                        for j=1:J
                            lock(lk) do
                                vJ[x,z,s,j,sex,t] = util_fn(x,z,s,j,sex, t) + xtran[x, :,z, j]' * FU[:, z, s, sex, t+1]
                            end
                        end
                        lock(lk) do
                            v0[x,z,s,sex,t] = xtran[x, :,z, 1]' * FU[:, z, s, sex, t+1]
                            FV[x,z, s, sex, t] =  (log(sum(exp.(vJ[x,z,s,:,sex,t])) + exp(v0[x,z,s,sex,t])) + eul)
                        end
                        for j=1:J
                            lock(lk) do
                                pJ_single[x,z,s,j,sex,t] = exp(vJ[x,z,s,j,sex,t])/(sum(exp.(vJ[x,z,s,:,sex,t]))+exp(v0[x,z,s,sex,t]))
                            end
                        end
                        lock(lk) do
                            p0_single[x,z,s,sex,t]=exp(v0[x,z,s,sex,t])/(sum(exp.(vJ[x,z,s,:,sex,t]))+exp(v0[x,z,s,sex,t]))
                        end
                        check=sum(pJ_single[x,z,s,:,sex,t])+p0_single[x,z,s,sex,t]
                        if abs(check-1)>1e-10
                            # println("CCPs don't sum up to 1! $check")
                        end
                    end
                end
            end
        end
        #############################
        # Choices as married
        #############################
        for sex in 1:2
            for s_m in 1:2
                for s_f in 1:2
                    for z_m in 1:zbin
                        for z_f in 1:zbin
                            for x_m in 1:xbin
                                Threads.@threads for x_f in 1:xbin
                                    for j_m=1:J+1
                                        for j_f=1:J+1
                                            # x_m1=Int(x_m*(j_m>1)+(j_m==1))
                                            # x_f1=Int(x_f*(j_f>1)+(j_f==1))
                                            # j_m1=Int((j_m-1)*(j_m>1)+(j_m==1))
                                            # j_f1=Int((j_f-1)*(j_f>1)+(j_f==1))
                                            lock(lk) do
                                                # v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,sex,t] = util_married_fn(x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,sex,t) + sum(sum(xtran_married[:,:, x_m1, x_f1, z_m, z_f, j_m1,j_f1] .* FU_married[:,:, z_m, z_f, s_m, s_f, sex, t+1],dims=(1,2)))
                                                v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,sex,t] = util_married_fn(x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,sex,t) + sum(sum(xtran_married[:,:, x_m, x_f, z_m, z_f, j_m,j_f] .* FU_married[:,:, z_m, z_f, s_m, s_f, sex, t+1],dims=(1,2)))
                                            end
                                        end
                                    end
                                    # Calculate CCP for couple's 2nd stage choices. Need to compute expected utilities
                                    denom=sum(exp.(lambda.*v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,1,t] .+ (1 .-lambda).*v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,2,t]))
                                    for j_m=1:J+1
                                        for j_f=1:J+1
                                            num=lambda.* v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,1,t] + (1 .- lambda) .* v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,2,t]
                                            lock(lk) do
                                                p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,j_m,j_f,t]=exp(num)/denom
                                            end
                                        end
                                    end
                                    check=sum(p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,t])
                                    if abs(check-1)>1e-15
                                        # println("CCPs don't sum up to 1! $check")
                                    end
                                    lock(lk) do
                                        FV_married[x_m,x_f,z_m,z_f,s_m,s_f,sex,t]=sum(p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,t].*(v0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,sex,t] .- log.(p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,t])))
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
            for s in 1:2
                for s_sp_orig in 1:2
                    for z in 1:zbin
                        for z_sp_orig in 1:zbin
                            Threads.@threads for x in 1:xbin 
                                for x_sp_orig in 1:xbin
                                    x_=xval[x]
                                    z_=zval[z]
                                    x_sp_orig_=xval[x_sp_orig]
                                    z_sp_orig_=zval[z_sp_orig]
                                    denom=0
                                    for s_sp_dest in 1:2
                                        for z_sp_dest in 1:zbin
                                            z_sp_dest_=zval[z_sp_dest]
                                            for x_sp_dest in 1:xbin
                                                x_sp_dest_=xval[x_sp_dest]
                                                lock(lk) do
                                                    u_marry_married[x,x_sp_orig,x_sp_dest,z,z_sp_orig,z_sp_dest,s,s_sp_orig,s_sp_dest,sex,t]=kappa(x_,x_sp_orig_,x_sp_dest_,z_,z_sp_orig_,z_sp_dest_,s,s_sp_orig,s_sp_dest,sex) + FV_married[x,x_sp_dest,z,z_sp_dest,s,s_sp_dest,sex,t]
                                                end
                                                denom=denom+exp(u_marry_married[x,x_sp_orig,x_sp_dest,z,z_sp_orig,z_sp_dest,s,s_sp_orig,s_sp_dest,sex,t] )
                                            end
                                        end
                                    end
                                    # Add the value of being single
                                    lock(lk) do
                                        u_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t]=kappa(x_,x_sp_orig_,0,z_,z_sp_orig_,0,s,s_sp_orig,0,sex) + FV[x,z, s, sex,t]
                                    end
                                    denom=denom+exp(u_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t])
                                    lock(lk) do
                                        FU_married[x,x_sp_orig, z, z_sp_orig, s, s_sp_orig, sex, t]= Beta* log(denom)
                                        p_marry_married[x,x_sp_orig,:,z,z_sp_orig,:,s,s_sp_orig,:,sex,t].=exp.(u_marry_married[x,x_sp_orig,:,z,z_sp_orig,:,s,s_sp_orig,:,sex,t])./denom
                                        p_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t]=exp(u_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t])/denom
                                    end
                                    check=sum(p_marry_married[x,x_sp_orig,:,z,z_sp_orig,:,s,s_sp_orig,:,sex,t])+p_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t]
                                    if abs(check-1)>1e-15
                                        # println("CCPs don't sum up to 1! $check")
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        # Expected value for people who are SINGLE at the beginning of stage 1
        for sex in 1:2
            for s in 1:2
                for z in 1:zbin
                    Threads.@threads for x in 1:xbin 
                        adj = x + (z - 1) * xbin
                        x_=xval[x]
                        z_=zval[z]
                        denom=0
                        for s_sp_dest in 1:2
                            for z_sp_dest in 1:zbin
                                z_sp_dest_=zval[z_sp_dest]
                                for x_sp_dest in 1:xbin
                                    x_sp_dest_=xval[x_sp_dest]
                                    lock(lk) do
                                        u_marry_single[x,x_sp_dest,z,z_sp_dest,s,s_sp_dest,sex,t]=kappa(x_,0,x_sp_dest_,z_,0,z_sp_dest_,s,0,s_sp_dest,sex) + FV_married[x,x_sp_dest,z,z_sp_dest,s,s_sp_dest,sex,t]
                                    end
                                    denom=denom+exp( u_marry_single[x,x_sp_dest,z,z_sp_dest,s,s_sp_dest,sex,t])
                                end
                            end
                        end
                        # Add the value of being single
                        lock(lk) do
                            u_single_single[x,z,s,sex,t]=kappa(x_,0,0,z_,0,0,s,0,0,sex) + FV[x,z, s, sex, t]
                        end
                        denom=denom + exp(u_single_single[x,z,s,sex,t])
                        lock(lk) do
                            FU[x, z, s, sex, t]= Beta* log(denom)
                            p_marry_single[x,:,z,:,s,:,sex,t].=exp.(u_marry_single[x,:,z,:,s,:,sex,t])./denom
                            p_single_single[x,z,s,sex,t]=exp(u_single_single[x,z,s,sex,t])/denom
                        end
                        check=sum(p_marry_single[x,:,z,:,s,:,sex,t])+p_single_single[x,z,s,sex,t]
                        if abs(check-1)>1e-15
                            # println("CCPs don't sum up to 1! $check")
                        end
                    end
                end
            end
        end
        t -= 1
        println(t)
    end
    toc=time()
    duration=toc-tic

    if true 
        # In SOEP, people enter the panel when they cohabit with someone who is already in the panel. I need to replicate this in the data generation
        # Each individual starts as single and I generate their history. At each period, there is an exogenous probability that an individual is sampled by SOEP, as a household.
        # If in that period the individual is married, then their spouse is now considered an individual too and I generate a history for them too.
        # After being sampled, at each period, there is an exogenous probability that and individual will drop out of the panel, for whatever reason.
        ID_1 = [(n, 0) for n in 1:N]
        # ID_1= collect(1:N)
        Sex_1 = (rand(N) .> .5) .+1
        # State is the unobserved type
        State_1 = Int.(rand(N) .> Pi) .+1
        State_sp_orig_1=Array{Union{Missing, Int64}}(undef, N, T2+1)
        State_sp_dest_1=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Y_1 = Array{Union{Missing, Int64}}(undef, N, T2)
        Y_sp_dest_1=Array{Union{Missing, Int64}}(undef, N, T2)
        Y_sp_orig_1=Array{Union{Missing, Int64}}(undef, N, T2)
        # Zstate is the randomly drawn position index on the grid of values of z
        Zstate_1 = Int.(ceil.(length(zval) * rand(N)))
        Zstate_sp_orig_1=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Zstate_sp_dest_1=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Z_1 = zval[Zstate_1]
        Z_sp_dest_1=Array{Union{Missing, Float64}}(undef, N, T2+1)
        Z_sp_orig_1=Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_1 = Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_1[:,1].=0
        X_sp_orig_1=Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_sp_dest_1=Array{Union{Missing, Float64}}(undef, N, T2+1)
        Xstate_1 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Xstate_1[:,1].=1
        Xstate_sp_orig_1 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Xstate_sp_dest_1 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Married_1=Array{Union{Missing, Int64}}(undef, N, T2)
        Sampled_1=zeros(Int,N,T2 )
        Draw_1_stage2 = rand(N, T2 )
        Draw_1_transition = rand(N, T2 )
        Draw_1_stage1 = rand(N, T2 )
        Draw_1_sample = rand(N, T2 )
        Draw_1_drop = rand(N, T2 )
        Temp_1=zeros(N,T2)

        ID_2=Array{Union{Missing, Tuple}}(missing, N)
        Sex_2 = Array{Union{Missing, Int64}}(undef, N)
        # State is the unobserved type
        State_2 = Array{Union{Missing, Int64}}(undef, N)
        State_sp_orig_2=Array{Union{Missing, Int64}}(undef, N, T2+1)
        State_sp_dest_2=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Y_2 = Array{Union{Missing, Int64}}(undef, N, T2)
        Y_sp_dest_2= Array{Union{Missing, Int64}}(undef, N, T2)
        Y_sp_orig_2= Array{Union{Missing, Int64}}(undef, N, T2)
        # Zstate is the randomly drawn position index on the grid of values of z
        Zstate_2 = Array{Union{Missing, Int64}}(undef, N)
        Zstate_sp_orig_2=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Zstate_sp_dest_2=Array{Union{Missing, Int64}}(undef, N, T2+1)
        Z_2 = Array{Union{Missing, Float64}}(undef, N)
        Z_sp_dest_2=Array{Union{Missing, Float64}}(undef, N, T2+1)
        Z_sp_orig_2=Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_2 = Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_sp_orig_2=Array{Union{Missing, Float64}}(undef, N, T2+1)
        X_sp_dest_2=Array{Union{Missing, Float64}}(undef, N, T2+1)
        Xstate_2 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Xstate_sp_orig_2 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Xstate_sp_dest_2 = Array{Union{Missing, Int64}}(undef, N, T2+1)
        Married_2=Array{Union{Missing, Int64}}(undef, N, T2)
        Sampled_2=zeros(Int,N,T2 )
        Draw_2_stage2 = rand(N, T2 )
        Draw_2_transition = rand(N, T2 )
        Draw_2_stage1 = rand(N, T2 )
        Draw_2_sample = rand(N, T2 )
        Draw_2_drop = rand(N, T2 )
        Temp_1=zeros(N,T2)

        temp1=reshape(collect(1:xbin*zbin*NTypes),(xbin,zbin,NTypes))
        temp2=reshape(collect(1:(J+1)^2),(J+1,J+1))
        temp3=reshape(collect(1:xbin_m*xbin_f),(xbin_m,xbin_f))

        for n in 1:N
            # Convert position index on Z space into Z value
            z=Zstate_1[n]
            z_ = zval[z]
            Z_1[n] = z_
            sex=Sex_1[n]
            # adj0 is the state achieved *immediately at t* if choose to replace engine at time t, i.e. it's zero mileage but ofc keep the same z value
            # But then, by t+1 you accumulate mileage
            # So, notice that this model is a bit nonstandard in that the choice also affects the current state. In other words, choice come before the realization of the current state.
            adj0 = (z - 1) * xbin + 1
            # z2 and z3 select the states that you can possibly achieve at t+1 given your value of z.
            # it's then up to xtran to assign a zero probability to states with lower mileage than your current mileage (which is reset to zero if you do replace your engine)
            z2 = (z - 1) * xbin + 1
            z3 = z2 + xbin - 1
            for t in 1:T2 
                # println(t)
                s=Int(State_1[n])
                x_=X_1[n, t]
                x=Xstate_1[n,t]
                adj = x + (z - 1) * xbin
                if t==1
                    married_tminus1=0
                else
                    married_tminus1=Married_1[n,t-1]
                    if married_tminus1==1
                        s_sp_orig=Int(State_sp_orig_1[n,t])
                        z_sp_orig=Zstate_sp_orig_1[n,t]
                        z_sp_orig_=Z_sp_orig_1[n,t]
                        x_sp_orig=Xstate_sp_orig_1[n,t]
                        x_sp_orig_=X_sp_orig_1[n,t]
                    end
                end

                ##########################
                # Stage 1: Choose partner
                ##########################
                if married_tminus1==1
                    p_marry=p_marry_married[x,x_sp_orig,:,z,z_sp_orig,:,s,s_sp_orig,:,sex,t]
                    p_single=p_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t]
                else
                    p_marry=p_marry_single[x,:,z,:,s,:,sex,t]
                    p_single=p_single_single[x,z,s,sex,t]
                end
                cumulprob=cumsum(vcat(p_single,p_marry[:]))
                spouse=searchsortedfirst(cumulprob,Draw_1_stage1[n, t])
                married=spouse>1

                Married_1[n,t]=married
                if married==1
                    # Because of how I did vcat above, spouse=1 means the person remains single
                    spouse=spouse-1
                    # Spouse now indicates the new spouse position in the array of size (xbin,zbin,NTypes)
                    Xstate_sp_dest_1[n,t]=findall(spouse.==temp1)[1][1]
                    X_sp_dest_1[n,t]=xval[Xstate_sp_dest_1[n,t]]
                    Zstate_sp_dest_1[n,t]=findall(spouse.==temp1)[1][2]
                    Z_sp_dest_1[n,t]=zval[Zstate_sp_dest_1[n,t]]
                    State_sp_dest_1[n,t]=findall(spouse.==temp1)[1][3]
                end

                ##########################
                # Stage 1.1: sampling
                ##########################
                P_sample=.3
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
                        x_m=x
                        x_f=Xstate_sp_dest_1[n,t]
                        z_m=z
                        z_f=Zstate_sp_dest_1[n,t]
                        s_m=s
                        s_f=State_sp_dest_1[n,t]
                    else
                        x_f=x
                        x_m=Xstate_sp_dest_1[n,t]
                        z_f=z
                        z_m=Zstate_sp_dest_1[n,t]
                        s_f=s
                        s_m=State_sp_dest_1[n,t]
                    end
                    p=p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,t] #Remember that in this array, j=1 means replacing the engine, it's the normalized action
                    temp=searchsortedfirst(cumsum(p[:]),Draw_1_stage2[n, t])
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
                    p=vcat(p0_single[x,z,s,sex,t],pJ_single[x,z,s,:,sex,t])
                    j=searchsortedfirst(cumsum(p[:]),Draw_1_stage2[n, t])
                    Y_1[n,t]=j
                end

                ###################################
                # State transition to next period
                ###################################
                if married==1
                    # j_m1=j_m*(j_m>0)+(j_m==0)
                    # x_m1=x_m*(j_m>0)+(j_m==0)
                    # j_f1=j_f*(j_f>0)+(j_f==0)
                    # x_f1=x_f*(j_f>0)+(j_f==0)
                    # temp=searchsortedfirst(xtranc_married[:,:,x_m1,x_f1,z_m,z_f,j_m1,j_f1][:],Draw_1_transition[n,t])
                    temp=searchsortedfirst(xtranc_married[:,:,x_m,x_f,z_m,z_f,j_m,j_f][:],Draw_1_transition[n,t])
                    Temp_1[n,t]=temp
                    # x_m=mod(temp,xbin)*(mod(temp,xbin)>0)+xbin*(mod(temp,xbin)==0)
                    # x_f=Int(ceil(temp/xbin))
                    position = findfirst(x -> x == temp, temp3)
                    x_m=position[1]
                    x_f=position[2]
                    if sex==1
                        Xstate_1[n, t+1] = x_m
                        X_1[n,t+1]=xval[x_m]
                        Xstate_sp_orig_1[n, t+1] = x_f
                        X_sp_orig_1[n,t+1]=x_f
                    else
                        Xstate_1[n, t+1] = x_f
                        X_1[n,t+1]=xval[x_f]
                        Xstate_sp_orig_1[n, t+1] = x_m
                        X_sp_orig_1[n,t+1]=xval[x_m]
                    end
                    State_sp_orig_1[n,t+1]=State_sp_dest_1[n,t]
                    Zstate_sp_orig_1[n,t+1]=Zstate_sp_dest_1[n,t]
                    Z_sp_orig_1[n,t+1]=Z_sp_dest_1[n,t]
                else
                    j_1=j*(j>0)+(j==0)
                    x_1=x*(j>0)+(j==0)
                    Xstate_1[n, t+1]=searchsortedfirst(xtranc[x_1,:,z,j_1],Draw_1_transition[n,t])
                    X_1[n, t+1] = xval[Xstate_1[n, t+1]]
                end
                P_drop=0.1
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
                # ID_2[n]=-ID_1[n]
                ID_2[n]=(-ID_1[n][1],0)
                # Convert position index on Z space into Z value
                z=Zstate_sp_dest_1[n,t_sample]
                z_ = zval[z]
                Z_2[n] = z_
                Sex_2[n]=(Sex_1[n]==1)*2+(Sex_1[n]==2)
                sex=Int(Sex_2[n])
                # adj0 is the state achieved *immediately at t* if choose to replace engine at time t, i.e. it's zero mileage but ofc keep the same z value
                # But then, by t+1 you accumulate mileage
                # So, notice that this model is a bit nonstandard in that the choice also affects the current state. In other words, choice come before the realization of the current state.
                adj0 = (z - 1) * xbin + 1
                # z2 and z3 select the states that you can possibly achieve at t+1 given your value of z.
                # it's then up to xtran to assign a zero probability to states with lower mileage than your current mileage (which is reset to zero if you do replace your engine)
                z2 = (z - 1) * xbin + 1
                z3 = z2 + xbin - 1
                s=State_sp_dest_1[n,t_sample]
                State_2[n]=s
                x_=X_sp_dest_1[n, t_sample]
                X_2[n,t_sample]=x_
                x=Xstate_sp_dest_1[n,t_sample]
                Xstate_2[n,t_sample]=x
                adj = x + (z - 1) * xbin
                Married_2[n,t_sample]=1
                married=1
                Xstate_sp_dest_2[n,t_sample]=Xstate_1[n,t_sample]
                X_sp_dest_2[n,t_sample]=xval[ Xstate_sp_dest_2[n,t_sample]]
                Zstate_sp_dest_2[n,t_sample]=Zstate_1[n]
                Z_sp_dest_2[n,t_sample]=zval[Zstate_sp_dest_2[n,t_sample]]
                State_sp_dest_2[n,t_sample]= State_1[n]
                for t in t_sample:T2 
                    Sampled_2[n,t]=1
                    if t>t_sample
                        x_=X_2[n, t]
                        x=Xstate_2[n,t]
                        adj = x + (z - 1) * xbin
                        married_tminus1=Married_2[n,t-1]
                        if married_tminus1==1
                            s_sp_orig=Int(State_sp_orig_2[n,t])
                            z_sp_orig=Zstate_sp_orig_2[n,t]
                            z_sp_orig_=Z_sp_orig_2[n,t]
                            x_sp_orig=Xstate_sp_orig_2[n,t]
                            x_sp_orig_=X_sp_orig_2[n,t]
                        end
                    
                        ##########################
                        # Stage 1: Choose partner
                        ##########################
                        if married_tminus1==1
                            p_marry=p_marry_married[x,x_sp_orig,:,z,z_sp_orig,:,s,s_sp_orig,:,sex,t]
                            p_single=p_single_married[x,x_sp_orig,z,z_sp_orig,s,s_sp_orig,sex,t]
                        else
                            p_marry=p_marry_single[x,:,z,:,s,:,sex,t]
                            p_single=p_single_single[x,z,s,sex,t]
                        end
                        cumulprob=cumsum(vcat(p_single,p_marry[:]))
                        spouse=searchsortedfirst(cumulprob,Draw_2_stage1[n, t])
                        married=spouse>1

                        Married_2[n,t]=married
                        if married==1
                            # Because of how I did vcat above, spouse=1 means the person remains single
                            spouse=spouse-1
                            # Spouse now indicates the new spouse position in the array of size (xbin,zbin,NTypes)
                            Xstate_sp_dest_2[n,t]=findall(spouse.==temp1)[1][1]
                            X_sp_dest_2[n,t]=xval[Xstate_sp_dest_2[n,t]]
                            Zstate_sp_dest_2[n,t]=findall(spouse.==temp1)[1][2]
                            Z_sp_dest_2[n,t]=zval[Zstate_sp_dest_2[n,t]]
                            State_sp_dest_2[n,t]=findall(spouse.==temp1)[1][3]
                        end

                    end

                    ##########################
                    # Stage 2: Other choices
                    ##########################
                    if married==1
                        if sex==1
                            x_m=x
                            x_f=Xstate_sp_dest_2[n,t]
                            z_m=z
                            z_f=Zstate_sp_dest_2[n,t]
                            s_m=s
                            s_f=State_sp_dest_2[n,t]
                        else
                            x_f=x
                            x_m=Xstate_sp_dest_2[n,t]
                            z_f=z
                            z_m=Zstate_sp_dest_2[n,t]
                            s_f=s
                            s_m=State_sp_dest_2[n,t]
                        end
                        p=p0J_married[x_m,x_f,z_m,z_f,s_m,s_f,:,:,t] #Remember that in this array, j=1 means replacing the engine, it's the normalized action
                        temp=searchsortedfirst(cumsum(p[:]),Draw_2_stage2[n, t])
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
                        p=vcat(p0_single[x,z,s,sex,t],pJ_single[x,z,s,:,sex,t])
                        j=searchsortedfirst(cumsum(p[:]),Draw_2_stage2[n, t])
                        Y_2[n,t]=j
                    end
                    

                    ###################################
                    # State transition to next period
                    ###################################
                    if married==1
                        # j_m1=j_m*(j_m>0)+(j_m==0)
                        # x_m1=x_m*(j_m>0)+(j_m==0)
                        # j_f1=j_f*(j_f>0)+(j_f==0)
                        # x_f1=x_f*(j_f>0)+(j_f==0)
                        # temp=searchsortedfirst(xtranc_married[:,:,x_m1,x_f1,z_m,z_f,j_m1,j_f1][:],Draw_2_transition[n,t])
                        temp=searchsortedfirst(xtranc_married[:,:,x_m,x_f,z_m,z_f,j_m,j_f][:],Draw_2_transition[n,t])
                        # x_m=mod(temp,xbin)*(mod(temp,xbin)>0)+xbin*(mod(temp,xbin)==0)
                        # x_f=Int(ceil(temp/xbin))
                        position = findfirst(x -> x == temp, temp3)
                        x_m=position[1]
                        x_f=position[2]
                        if sex==1
                            Xstate_2[n, t+1] = x_m
                            X_2[n,t+1]=xval[x_m]
                            Xstate_sp_orig_2[n, t+1] = x_f
                            X_sp_orig_2[n,t+1]=x_f
                        else
                            Xstate_2[n, t+1] = x_f
                            X_2[n,t+1]=xval[x_f]
                            Xstate_sp_orig_2[n, t+1] = x_m
                            X_sp_orig_2[n,t+1]=xval[x_m]
                        end
                        State_sp_orig_2[n,t+1]=State_sp_dest_2[n,t]
                        Zstate_sp_orig_2[n,t+1]=Zstate_sp_dest_2[n,t]
                        Z_sp_orig_2[n,t+1]=Z_sp_dest_2[n,t]
                    else
                        j_2=j*(j>0)+(j==0)
                        x_2=x*(j>0)+(j==0)
                        Xstate_2[n, t+1]=searchsortedfirst(xtranc[x_2,:,z,j_2],Draw_2_transition[n,t])
                        X_2[n, t+1] = xval[Xstate_2[n, t+1]]
                    end

                    P_drop=0.05
                    Dropped=rand()<P_drop
                    if Dropped==1
                        break
                    end

                end
            end
        end


        sampled_1 = findall(row -> any(x -> x == 1, row), eachrow(Sampled_1))
        sampled_2 = maximum(Sampled_2[sampled_1,:], dims=2)
        

        Sampled_1=Sampled_1[sampled_1,1:T2]
        Sampled_2=Sampled_2[sampled_1,1:T2]
        FirstSampled_1= [findfirst(x -> x == 1, row) for row in eachrow(Sampled_1)][:]
        Y_1=ifelse.(Sampled_1.==0,missing,Y_1[sampled_1,1:T2])
        Y_sp_dest_1=ifelse.(Sampled_1.==0,missing,Y_sp_dest_1[sampled_1,1:T2])
        Y_2=ifelse.(Sampled_2.==0,missing,Y_2[sampled_1,1:T2])
        Y_sp_dest_2=ifelse.(Sampled_2.==0,missing,Y_sp_dest_2[sampled_1,1:T2])
        X_1=ifelse.(Sampled_1.==0,missing,X_1[sampled_1,1:T2])
        X_sp_dest_1=ifelse.(Sampled_1.==0,missing,X_sp_dest_1[sampled_1,1:T2])
        X_sp_orig_1=ifelse.(Sampled_1.==0,missing,X_sp_orig_1[sampled_1,1:T2])
        X_2=ifelse.(Sampled_2.==0,missing,X_2[sampled_1,1:T2])
        X_sp_dest_2=ifelse.(Sampled_2.==0,missing,X_sp_dest_2[sampled_1,1:T2])
        X_sp_orig_2=ifelse.(Sampled_2.==0,missing,X_sp_orig_2[sampled_1,1:T2])
        Xstate_1=ifelse.(Sampled_1.==0,missing,Xstate_1[sampled_1,1:T2])
        Xstate_sp_dest_1=ifelse.(Sampled_1.==0,missing,Xstate_sp_dest_1[sampled_1,1:T2])
        Xstate_sp_orig_1=ifelse.(Sampled_1.==0,missing,Xstate_sp_orig_1[sampled_1,1:T2])
        Xstate_2=ifelse.(Sampled_2.==0,missing,Xstate_2[sampled_1,1:T2])
        Xstate_sp_dest_2=ifelse.(Sampled_2.==0,missing,Xstate_sp_dest_2[sampled_1,1:T2])
        Xstate_sp_orig_2=ifelse.(Sampled_2.==0,missing,Xstate_sp_orig_2[sampled_1,1:T2])
        Z_1=Z_1[sampled_1,:]
        Z_sp_dest_1=ifelse.(Sampled_1.==0,missing,Z_sp_dest_1[sampled_1,1:T2])
        Z_sp_orig_1=ifelse.(Sampled_1.==0,missing,Z_sp_orig_1[sampled_1,1:T2])
        Z_2=ifelse.(sampled_2.==0,missing,Z_2[sampled_1,:])
        Z_sp_dest_2=ifelse.(Sampled_2.==0,missing,Z_sp_dest_2[sampled_1,1:T2])
        Z_sp_orig_2=ifelse.(Sampled_2.==0,missing,Z_sp_orig_2[sampled_1,1:T2])
        Zstate_1=Zstate_1[sampled_1,:]
        Zstate_sp_dest_1=ifelse.(Sampled_1.==0,missing,Zstate_sp_dest_1[sampled_1,1:T2])
        Zstate_sp_orig_1=ifelse.(Sampled_1.==0,missing,Zstate_sp_orig_1[sampled_1,1:T2])
        Zstate_2=ifelse.(sampled_2.==0,missing,Zstate_2[sampled_1,:])
        Zstate_sp_dest_2=ifelse.(Sampled_2.==0,missing,Zstate_sp_dest_2[sampled_1,1:T2])
        Zstate_sp_orig_2=ifelse.(Sampled_2.==0,missing,Zstate_sp_orig_2[sampled_1,1:T2])
        Sex_1=Sex_1[sampled_1,:]
        Sex_2=ifelse.(sampled_2.==0,missing,Sex_2[sampled_1,:])
        ID_1=ID_1[sampled_1,:]
        ID_2=ifelse.(sampled_2.==0,missing,ID_2[sampled_1])
        Married_1=ifelse.(Sampled_1.==0,missing,Married_1[sampled_1,1:T2])
        Married_2=ifelse.(Sampled_2.==0,missing,Married_2[sampled_1,1:T2])

        N=size(Sampled_1)[1]

        # Additionally, all _sp_orig variables should have missing values in the first period an agent is in the sample
        for n in 1:N
            tfirst=FirstSampled_1[n]
            X_sp_orig_1[n,tfirst]=missing
            X_sp_orig_2[n,tfirst]=missing
            Xstate_sp_orig_1[n,tfirst]=missing
            Xstate_sp_orig_2[n,tfirst]=missing
            Z_sp_orig_1[n,tfirst]=missing
            Z_sp_orig_2[n,tfirst]=missing
            Zstate_sp_orig_1[n,tfirst]=missing
            Zstate_sp_orig_2[n,tfirst]=missing
        end
    end


    # Create ID for spouses too and mapping from individuals to spouse for each t
    ID_sp_dest_1=Array{Union{Missing, Tuple}}(missing, N, T2)
    ID_sp_dest_2=Array{Union{Missing, Tuple}}(missing, N, T2)
    
    for n in 1:N
        id_sp_1=ifelse(ismissing(ID_2[n]),ID_1[n],ID_2[n])
        id_sp_2=ID_1[n]
        condition=false
        for t in 1:T2
            if t<FirstSampled_1[n]
                continue
            end
            if (t==FirstSampled_1[n]) & (Married_1[n,t]===1)
                ID_sp_dest_1[n,t]=id_sp_1
                ID_sp_dest_2[n,t]=id_sp_2
            else
                condition_1=(Xstate_sp_orig_1[n,t]!==Xstate_sp_dest_1[n,t]) | (Zstate_sp_orig_1[n,t]!==Zstate_sp_dest_1[n,t]) | (State_sp_orig_1[n,t]!==State_sp_dest_1[n,t])
                condition_2=(Xstate_sp_orig_2[n,t]!==Xstate_sp_dest_2[n,t]) | (Zstate_sp_orig_2[n,t]!==Zstate_sp_dest_2[n,t]) | (State_sp_orig_2[n,t]!==State_sp_dest_2[n,t])
                if t>1
                    condition=(ID_sp_dest_1[n,t-1]===ID_2[n]) & (ID_sp_dest_1[n,t-1]!==missing) & (condition_1 | condition_2==true)
                end
                if condition==true
                    id_sp_1=(ID_1[n][1],0)
                    id_sp_2=(ID_2[n][1],0)
                end
                if (condition_1==true) & (Married_1[n,t]===1)
                    id_sp_1=(id_sp_1[1],id_sp_1[2]+1)
                end
                if (condition_2==true) & (Married_2[n,t]===1)
                    id_sp_2=(id_sp_2[1],id_sp_2[2]+1)
                end
                if (Married_1[n,t]===1)
                    ID_sp_dest_1[n,t]=id_sp_1
                end
                if (Married_2[n,t]===1)
                    ID_sp_dest_2[n,t]=id_sp_2
                end
            end
            
        end
    end



    dt_mat=Data_Matrices(Y_1, 
    Y_sp_dest_1, 
    Y_sp_orig_1,
    X_1, 
    X_sp_dest_1, 
    X_sp_orig_1, 
    Z_1, 
    Z_sp_dest_1,
    Z_sp_orig_1,
    Xstate_1, 
    Xstate_sp_dest_1, 
    Xstate_sp_orig_1, 
    Zstate_1, 
    Zstate_sp_dest_1, 
    Zstate_sp_orig_1, 
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
    X_2, 
    X_sp_dest_2, 
    X_sp_orig_2, 
    Z_2, 
    Z_sp_dest_2, 
    Z_sp_orig_2,
    Xstate_2, 
    Xstate_sp_dest_2, 
    Xstate_sp_orig_2, 
    Zstate_2, 
    Zstate_sp_dest_2, 
    Zstate_sp_orig_2, 
    State_2, 
    State_sp_dest_2, 
    State_sp_orig_2, 
    ID_2, 
    Sampled_2, 
    Married_2, 
    Sex_2, 
    Adj,
    ID_sp_dest_1,
    ID_sp_dest_1,
    FirstSampled_1)

    return dt_mat
end

struct Data_Matrices
    Y_1:: Array 
    Y_sp_dest_1:: Array 
    Y_sp_orig_1:: Array 
    X_1:: Array 
    X_sp_dest_1:: Array 
    X_sp_orig_1:: Array 
    Z_1:: Array 
    Z_sp_dest_1:: Array
    Z_sp_orig_1:: Array
    Xstate_1:: Array 
    Xstate_sp_dest_1:: Array 
    Xstate_sp_orig_1:: Array 
    Zstate_1:: Array 
    Zstate_sp_dest_1:: Array 
    Zstate_sp_orig_1:: Array 
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
    X_2:: Array 
    X_sp_dest_2:: Array 
    X_sp_orig_2:: Array 
    Z_2:: Array 
    Z_sp_dest_2:: Array 
    Z_sp_orig_2:: Array
    Xstate_2:: Array 
    Xstate_sp_dest_2:: Array 
    Xstate_sp_orig_2:: Array 
    Zstate_2:: Array 
    Zstate_sp_dest_2:: Array 
    Zstate_sp_orig_2:: Array 
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
