function generate_data_wrap(load_=false)
    
    if load_==false
        dm = genbus4(alpha, N_, T, xtran, xtranc, xtran_married, xtranc_married, xbin, zbin, xval, zval, T2)
        # save("busdata0210_dataset.jld", "dm", dm, "xtran", xtran, "xtranc", xtranc, "xtran_married", xtran_married, 
        # "N_", N_, "xval", xval, "zval", zval, "xbin", xbin, "zbin", zbin, "T2", T2)
    else
        vars = load("busdata0210_dataset.jld")
        for (key, value) in vars
            @eval $(Symbol(key)) = $value
            eval(Meta.parse("$(key) = $(value)"))
        end
        # dm.Zstate_1=Int.(dm.Zstate_1)
        # dm.Xstate_1=Int.(dm.Xstate_1)
        # dm.Y_1=Int.(dm.Y_1)
    end

    if N_<length(dm.ID_1)
        # If the loaded dataset has many individuals and want to play with a smaller sample N_
        dm.ID_1=dm.ID_1[1:N_]
        dm.Zstate_1=dm.Zstate_1[1:N_]
        dm.Zstate_sp_dest_1=dm.Zstate_sp_dest_1[1:N_,:]
        dm.Zstate_sp_orig_1=dm.Zstate_sp_orig_1[1:N_,:]
        dm.Xstate_1=dm.Xstate_1[1:N_,:]
        dm.Xstate_sp_dest_1=dm.Xstate_sp_dest_1[1:N_,:]
        dm.Xstate_sp_orig_1=dm.Xstate_sp_orig_1[1:N_,:]
        dm.State_sp_dest_1=dm.State_sp_dest_1[1:N_,:]
        dm.State_sp_orig_1=dm.State_sp_orig_1[1:N_,:]
        dm.Y_1=dm.Y_1[1:N_,:]
        dm.Y_sp_dest_1=dm.Y_sp_dest_1[1:N_,:]
        dm.Y_sp_orig_1=dm.Y_sp_orig_1[1:N_,:]
        dm.X_1=dm.X_1[1:N_,:]
        dm.X_sp_dest_1=dm.X_sp_dest_1[1:N_,:]
        dm.X_sp_orig_1=dm.X_sp_orig_1[1:N_,:]
        dm.Z_1=dm.Z_1[1:N_]
        dm.Z_sp_dest_1=dm.Z_sp_dest_1[1:N_]
        dm.Z_sp_orig_1=dm.Z_sp_orig_1[1:N_]
        dm.Sampled_1=dm.Sampled_1[1:N_,:]
        dm.Married_1=dm.Married_1[1:N_,:]
        dm.Sex_1=dm.Sex_1[1:N_]
        dm.ID_2=dm.ID_2[1:N_]
        dm.Zstate_2=dm.Zstate_2[1:N_]
        dm.Zstate_sp_dest_2=dm.Zstate_sp_dest_2[1:N_,:]
        dm.Zstate_sp_orig_2=dm.Zstate_sp_orig_2[1:N_,:]
        dm.Xstate_2=dm.Xstate_2[1:N_,:]
        dm.X_sp_dest_2=dm.X_sp_dest_2[1:N_,:]
        dm.Xstate_sp_orig_2=dm.Xstate_sp_orig_2[1:N_,:]
        dm.State_sp_dest_2=dm.State_sp_dest_2[1:N_,:]
        dm.State_sp_orig_2=dm.State_sp_orig_2[1:N_,:]
        dm.Y_2=dm.Y_2[1:N_,:]
        dm.Y_sp_dest_2=dm.Y_sp_dest_2[1:N_,:]
        dm.Y_sp_orig_2=dm.Y_sp_orig_2[1:N_,:]
        dm.X_2=dm.X_2[1:N_,:]
        dm.X_sp_dest_2=dm.X_sp_dest_2[1:N_,:]
        dm.X_sp_orig_2=dm.X_sp_orig_2[1:N_,:]
        dm.Z_2=dm.Z_2[1:N_]
        dm.Z_sp_dest_2=dm.Z_sp_dest_2[1:N_]
        dm.Z_sp_orig_2=dm.Z_sp_orig_2[1:N_]
        dm.Sampled_2=dm.Sampled_2[1:N_,:]
        dm.Married_2=dm.Married_2[1:N_,:]
        dm.Sex_2=dm.Sex_2[1:N_]
        dm.State_1=dm.State_1[1:N_]
        dm.State_2=dm.State_2[1:N_]
    end

    N=length(dm.ID_1)

    ###########################################################
    # !!!!!!!!!! MUST READ THIS !!!!!!!!!!!!
    ###########################################################
    # In the following sections I create the vectors of variables used to create the vectors of residuals δ and instruments X used to  generate sample analogues of the moment conditions E(Xδ) for GMM estimation. The very same vectors are used to construct likelihoods and pseudo-likelihoods as well.
    # Notice that I form the vectors of observations so that they are grouped together for each period. The first 2N entries in any vector contain observations from t=1. The first N observations pertain to agent _1, the following N observations pertain to agent _2. Then, I stack these observations across time, i.e. the entries from 2N+1 to 2*2N pertain to t=2. And so on.
    # Then, I duplicate these observations as many times as there are combinations of each agent's unobserved type, and stack them vertically.
    # ALL VECTORS THAT CONTAIN OBSERVED DATA IN THIS ENTIRE MONTE CARLO EXERCISE ARE BUILT FOLLOWING THIS ORDERING
    # That each vector of residuals follow the same order is necessary for estimation because there might be (depending the moment considered) correlation between X_{it}δ_{it}(a) and X_{it}δ_{it}(a') for any individual i at time t and any two choices a and a', or between X_{ijt}δ_{ijt}(a) X_{ijt}δ_{ijt}(a') for any couple ij at time t and any two choices a and a', and that correlation needs to be taken into account when estimating standard errors.

    ###############################################
    # Spouses identifiers identifiers
    ###############################################
    # Notice that we always know what individual each observation pertains to based on their position in the vector
    id_sp=vcat(dm.ID_sp_dest_1,dm.ID_sp_dest_2)

    ###########################################################
    # Indicator for when a person is in the sample or not
    ###########################################################
    # Each originally sampled agent faces an exogenous sample drop-out probability
    # Furthermore, spouses of originally sample agents drop out of the sample if they are not married to them anymore.

    Insample=vcat(dm.Sampled_1,dm.Sampled_2)
    insample=Insample[:]
    insample=repeat(insample,NTypes^2)

    ###########################################################
    #   Indicator for first observation of an individual
    ###########################################################
    firstsampled=zeros(Int,size(Insample))
    for n in 1:N
        col=dm.FirstSampled_1[n]
        firstsampled[n,col]=1
        firstsampled[N+n,col]=1
    end
    firstsampled=firstsampled[:]
    firstsampled=repeat(firstsampled,NTypes^2)


    ###########################################################
    # Indicator for when a person is married or not
    ###########################################################
    Married=vcat(dm.Married_1,dm.Married_2)
    married=sparse(repeat(coalesce.(Married[:], NaN), NTypes^2)) # sparse(repeat(Married[:],NTypes^2))

    Married_tminus1_1=hcat(Array{Union{Missing, Int64}}(missing, N),dm.Married_1[:,2:end])
    Married_tminus1_2=hcat(Array{Union{Missing, Int64}}(missing, N),dm.Married_2[:,2:end])
    Married_tminus1=vcat(Married_tminus1_1,Married_tminus1_2)
    married_tminus1=sparse(repeat(coalesce.(Married_tminus1[:], NaN), NTypes^2)) #sparse(repeat(Married_tminus1[:],NTypes^2))

    ###########################################################
    # Sex
    ###########################################################
    sex=repeat(vcat(dm.Sex_1,dm.Sex_2),T2)[:]
    sex=repeat(sex,NTypes^2)

    ###########################################################
    # Spouse id
    ###########################################################
    id=vcat(dm.ID_sp_dest_1,dm.ID_sp_dest_2)[:]

    ###########################################################
    # Regressors for transitions estimation
    ###########################################################
    
    # Agents' transitions
    X_dest=vcat(dm.X_1[:,2:end],dm.X_2[:,2:end])
    X_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_dest)
    X_orig=vcat(dm.X_1[:,1:end-1],dm.X_2[:,1:end-1])
    X_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_orig)
    Xstate_dest=vcat(dm.Xstate_1[:,2:end],dm.Xstate_2[:,2:end])
    Xstate_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Xstate_dest)
    Xstate_orig=vcat(dm.Xstate_1[:,1:end-1],dm.Xstate_2[:,1:end-1])
    Xstate_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Xstate_orig)
    
    Z_orig=vcat(repeat(dm.Z_1,1,T2)[:,2:end],repeat(dm.Z_2,1,T2)[:,2:end])
    Z_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Z_orig)
    Zstate_orig=vcat(repeat(dm.Zstate_1,1,T2)[:,2:end],repeat(dm.Zstate_2,1,T2)[:,2:end])
    Zstate_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Zstate_orig)

    Y_orig=vcat(dm.Y_1[:,1:end-1],dm.Y_2[:,1:end-1])
    Y_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Y_orig)
    
    x_orig_fhat=repeat(X_orig[:],NTypes^2)
    x_dest_fhat=repeat(X_dest[:],NTypes^2)
    xstate_orig_fhat=repeat(Xstate_orig[:],NTypes^2)
    xstate_dest_fhat=repeat(Xstate_dest[:],NTypes^2)
    z_orig_fhat=repeat(Z_orig[:],NTypes^2)
    zstate_orig_fhat=repeat(Zstate_orig[:],NTypes^2)
    y_orig_fhat=repeat(Y_orig[:],NTypes^2)

    # Remember that by construction, replacing the engine (y=1) instantaneously sets mileage to zero, so
    x_orig_fhat[y_orig_fhat.===1].=0
    xstate_orig_fhat[y_orig_fhat.===1].=1

    # Spouses' transitions
    # The matrix X_sp_dest is built using X_sp_orig_1 because in X_sp_orig_1 "orig" refers to the origin spouse, i.e. the spouse chosen at t-1, but when thinking of state transitions, X_sp_orig_1 contains the destination *state* of the spouse as they transition from t-1 to t.
    X_sp_dest=vcat(dm.X_sp_orig_1[:,2:end],dm.X_sp_orig_2[:,2:end])
    X_sp_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_sp_dest)
    X_sp_orig=vcat(dm.X_sp_dest_1[:,1:end-1],dm.X_sp_dest_2[:,1:end-1])
    X_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_sp_orig)
    Xstate_sp_dest=vcat(dm.Xstate_sp_orig_1[:,2:end],dm.Xstate_sp_orig_2[:,2:end])
    Xstate_sp_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Xstate_sp_dest)
    Xstate_sp_orig=vcat(dm.Xstate_sp_dest_1[:,1:end-1],dm.Xstate_sp_dest_2[:,1:end-1])
    Xstate_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Xstate_sp_orig)
    
    Z_sp_orig=vcat(dm.Z_sp_orig_1[:,2:end],dm.Z_sp_orig_2[:,2:end])
    Z_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Z_sp_orig)
    Zstate_sp_orig=vcat(dm.Zstate_sp_orig_1[:,2:end],dm.Zstate_sp_orig_2[:,2:end])
    Zstate_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Zstate_sp_orig)

    Y_sp_orig=vcat(dm.Y_sp_dest_1[:,1:end-1],dm.Y_sp_dest_2[:,1:end-1])
    Y_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Y_sp_orig)
    
    x_sp_orig_fhat=repeat(X_sp_orig[:],NTypes^2)
    x_sp_dest_fhat=repeat(X_sp_dest[:],NTypes^2)
    xstate_sp_orig_fhat=repeat(Xstate_sp_orig[:],NTypes^2)
    xstate_sp_dest_fhat=repeat(Xstate_sp_dest[:],NTypes^2)
    z_sp_orig_fhat=repeat(Z_sp_orig[:],NTypes^2)
    zstate_sp_orig_fhat=repeat(Zstate_sp_orig[:],NTypes^2)
    y_sp_orig_fhat=repeat(Y_sp_orig[:],NTypes^2)

    # Remember that by construction, replacing the engine (y=1) instantaneously sets mileage to zero, so
    x_sp_orig_fhat[y_sp_orig_fhat.===1].=0
    xstate_sp_orig_fhat[y_sp_orig_fhat.===1].=1

    ###########################################################
    # Regressors for second-stage choice probabilities
    ###########################################################

    # Individual characteristics for choice probabilities as single
    # spouses by definition are missing from the sample if agent 1 or 2 are single.
    # These very same vectors can be used as regressors for first-stage choice probabilities, since the decision maker is agent 1 or 2 and not their spouses, so it is appropriate that there be missing values in correspondence of sp_1 and sp_2
    x=vcat(dm.X_1,dm.X_2)[:] # The reason I just repeat X_1 and X_2 is that only agent 1 and 2's choices enter the likelihood, not their spouses. Repeating X_! and X_@ allows calculating covariances with transitions 
    xstate=vcat(dm.Xstate_1,dm.Xstate_2)[:]
    Z_1_temp=ifelse.(dm.Sampled_1.==0,missing,repeat(dm.Z_1,1,T2))
    Z_2_temp=ifelse.(dm.Sampled_2.==0,missing,repeat(dm.Z_2,1,T2))
    z=vcat(Z_1_temp,Z_2_temp)[:]
    Zstate_1_temp=ifelse.(dm.Sampled_1.==0,missing,repeat(dm.Zstate_1,1,T2))
    Zstate_2_temp=ifelse.(dm.Sampled_2.==0,missing,repeat(dm.Zstate_2,1,T2))
    zstate=vcat(Zstate_1_temp,Zstate_2_temp)[:]
    y=vcat(dm.Y_1,dm.Y_2)[:]

    x=repeat(x,NTypes^2)
    xstate=repeat(xstate,NTypes^2)
    z=repeat(z,NTypes^2)
    zstate=repeat(zstate,NTypes^2)
    y=repeat(y,NTypes^2)


    # Couple characteristics for choice probabilities as married
    ## X
    x_sp_temp=vcat(dm.X_sp_dest_1,dm.X_sp_dest_2)[:]
    x_sp_temp=repeat(x_sp_temp,NTypes^2)
    xstate_sp_temp=vcat(dm.Xstate_sp_dest_1,dm.Xstate_sp_dest_2)[:]
    xstate_sp_temp=repeat(xstate_sp_temp,NTypes^2)

    x_m_cpl= x .* (sex.==1) .+ x_sp_temp .* (sex.==2)
    x_f_cpl= x .* (sex.==2) .+ x_sp_temp .* (sex.==1)

    xstate_m_cpl= xstate .* (sex.==1) + xstate_sp_temp .* (sex.==2)
    xstate_f_cpl= xstate .* (sex.==2) + xstate_sp_temp .* (sex.==1)

    ## Z
    z_sp_temp=vcat(dm.Z_sp_dest_1,dm.Z_sp_dest_2)[:]
    z_sp_temp=repeat(z_sp_temp,NTypes^2)
    zstate_sp_temp=vcat(dm.Zstate_sp_dest_1,dm.Zstate_sp_dest_2)[:]
    zstate_sp_temp=repeat(zstate_sp_temp,NTypes^2)

    z_m_cpl= z .* (sex.==1) .+ z_sp_temp .* (sex.==2)
    z_f_cpl= z .* (sex.==2) .+ z_sp_temp .* (sex.==1)

    zstate_m_cpl= zstate .* (sex.==1) + zstate_sp_temp .* (sex.==2)
    zstate_f_cpl= zstate .* (sex.==2) + zstate_sp_temp .* (sex.==1)

    ## Y
    y_sp_temp=vcat(dm.Y_sp_dest_1,dm.Y_sp_dest_2)[:]
    y_sp_temp=repeat(y_sp_temp,NTypes^2)

    y_m_cpl= y .* (sex.==1) .+ y_sp_temp .* (sex.==2)
    y_f_cpl= y .* (sex.==2) .+ y_sp_temp .* (sex.==1)

    ###########################################################
    # Regressors for first-stage choice probabilities
    ###########################################################

    x_sp_orig=vcat(dm.X_sp_orig_1,dm.X_sp_orig_2)[:]
    x_sp_orig=repeat(x_sp_orig,NTypes^2)
    xstate_sp_orig=vcat(dm.Xstate_sp_orig_1,dm.Xstate_sp_orig_2)[:]
    xstate_sp_orig=repeat(xstate_sp_orig,NTypes^2)
    x_sp_dest=vcat(dm.X_sp_dest_1,dm.X_sp_dest_2)[:]
    x_sp_dest=repeat(x_sp_dest,NTypes^2)
    xstate_sp_dest=vcat(dm.Xstate_sp_dest_1,dm.Xstate_sp_dest_2)[:]
    xstate_sp_dest=repeat(xstate_sp_dest,NTypes^2)

    z_sp_orig=vcat(dm.Z_sp_orig_1,dm.Z_sp_orig_2)[:]
    z_sp_orig=repeat(z_sp_orig,NTypes^2)
    zstate_sp_orig=vcat(dm.Zstate_sp_orig_1,dm.Zstate_sp_orig_2)[:]
    zstate_sp_orig=repeat(zstate_sp_orig,NTypes^2)
    z_sp_dest=vcat(dm.Z_sp_dest_1,dm.Z_sp_dest_2)[:]
    z_sp_dest=repeat(z_sp_dest,NTypes^2)
    zstate_sp_dest=vcat(dm.Zstate_sp_dest_1,dm.Zstate_sp_dest_2)[:]
    zstate_sp_dest=repeat(zstate_sp_dest,NTypes^2)

    ###############################################################
    # Placeholder vector of the unobserved state variable
    ###############################################################
    # Here's where we generate the "fake" matrix of regressors s. Remember that for each period the first N observations are individual 1, the second N are individual 1's spouse, the next N observations are individual 2 and the last N observations are for individual 2's spouse. To create all possible combinations, I count from 1 to 16 in binary with the following matrix. Then, I use kron to create the vector of regressors
    # I also need an indicator that only records the type of the destination spouse and the origin spouse
    s_=zeros(Int,1,NTypes^2)
    s_sp_orig_=zeros(Int,1,NTypes^2)
    s_sp_dest_=zeros(Int,1,NTypes^2)
    col=1
    # for j_sp_dest=1:NTypes
        # for j_sp_orig=1:NTypes
            # for j=1:NTypes
                for i_sp_dest=1:NTypes
                    # for i_sp_orig=1:NTypes
                        for i=1:NTypes
                            s_[:,col]=[i]
                            # s_sp_orig_[:,col]=[i_sp_orig,i,j_sp_orig,j]
                            s_sp_dest_[:,col]=[i_sp_dest]
                            col=col+1
                        end
                    # end
                end
            # end
        # end
    # end
    s_t=kron(s_,ones(Int,2*N))
    # s_t_sp_orig=kron(s_sp_orig_,ones(Int,N))
    s_t_sp_dest=kron(s_sp_dest_,ones(Int,2*N))
    s = repeat(s_t , T2,1)[:]


    ############################
    # Used in first step CCPs
    ############################
    s_sp_orig=ones(Int, size(s))
    s_sp_dest=repeat(s_t_sp_dest , T2,1)[:]
    
    ############################
    # Used in second step CCPs
    ############################
    s_m_cpl = s.*(sex.===1) .+ s_sp_dest.*(sex.===2)
    s_m_cpl = ifelse.(((married.===0) .| (insample.==0))[:], NaN, s_m_cpl)
    s_f_cpl = s.*(sex.===2) .+ s_sp_dest.*(sex.===1)
    s_f_cpl = ifelse.(((married.===0) .| (insample.==0))[:], NaN, s_f_cpl)
    
    ########################################################################
    # Time fixed effects
    ########################################################################
    # t2 simply indicates how far each time t is in the length of the panel
    t = kron(collect(1:T2), ones(Int, 2*N))
    t = repeat(t, NTypes^2)
    t2=t./T2
    
    # td are time period dummies to associate to the raw data
    td = zeros(NTypes^2 * 2*N * T2, T2 - 1)
    t_ = 1
    while t_ < T2
        td[:, t_] = (t2 .* T2 .- 1) .== t_
        t_ += 1
    end
    td=sparse(td)


    ##########################################################################
    #   Choice indicators, used as outcome variables in GMM estimators
    ##########################################################################
    # For every individual, generate a set of indicator variables for the first-stage and second-stage choices
    # First-stage: one indicator for each type of partner agents can choose
    I_st1=Array{Union{Missing, Float64}}(missing, size(x)[1],xbin_sp_dest,zbin_sp_dest,NTypes_sp_dest)
    for obs in 1:size(x)[1]
        for xstate_sp in 1:xbin_sp_dest
            for zstate_sp in 1:zbin_sp_dest
                for s_sp in 1:NTypes_sp_dest
                    # Triple inequality because these vectors might contain missing values
                    if ismissing(xstate_sp_dest[obs].*zstate_sp_dest[obs].*s_sp_dest[obs])
                        continue
                    end
                    I_st1[obs,xstate_sp,zstate_sp,s_sp]=(xstate_sp_dest[obs]===xstate_sp_dest)*(zstate_sp_dest[obs]===zstate_sp_dest)*(s_sp_dest[obs]===s_sp_dest)
                end
            end
        end
    end
    # I_st1 will need to be comformable to element-wise subtraction with the CCP arrays CCP_struct.p_st2_f and CCP_struct.p_st2_m, so I need to add some placeholder 

    # Second-stage: one set of indicators for each choice the agent can choose as single, and one set of indicators for each choice a couple can choose
    I_st2_single=Array{Union{Missing, Float64}}(missing, size(x)[1],J+1)
    for obs in 1:size(x)[1]
        for a in 1:J+1
            # Triple inequality because these vectors might contain missing values
            if ismissing(y[obs])
                continue
            end
            I_st2_single[obs,a]=(y[obs]===a)
        end
    end

    I_st2_cpl=Array{Union{Missing, Float64}}(missing, size(x)[1],J_m+1,J_f+1)
    for obs in 1:size(x)[1]
        for a_m in 1:J+1
            for a_f in 1:J+1
                # Triple inequality because these vectors might contain missing values
                if ismissing(y_m_cpl[obs]*y_f_cpl[obs])
                    continue
                end
                I_st2_cpl[obs,a_m,a_f]=(y_m_cpl[obs]===a_m).*(y_f_cpl[obs]===a_f)
            end
        end
    end


  

    data=Data(dm.Sex_1,
    dm.X_1,
    dm.Z_1,
    dm.State_1,
    dm.Zstate_1,
    dm.Xstate_1,
    dm.Y_1,
    dm.Married_1,
    dm.Sampled_1,
    dm.Sex_2,
    dm.X_2, 
    dm.Z_2, 
    dm.State_2, 
    dm.Zstate_2, 
    dm.Xstate_2, 
    dm.Y_2,
    dm.Married_2,
    dm.Sampled_2,
    dm.Zstate_sp_orig_1, 
    dm.Zstate_sp_dest_1, 
    dm.Xstate_sp_orig_1, 
    dm.Xstate_sp_dest_1, 
    dm.Y_sp_dest_1,
    dm.Zstate_sp_orig_2, 
    dm.Zstate_sp_dest_2, 
    dm.Xstate_sp_orig_2, 
    dm.Xstate_sp_dest_2, 
    dm.Y_sp_dest_2,
    dm.X_sp_orig_1, 
    dm.X_sp_dest_1, 
    dm.Z_sp_orig_1, 
    dm.Z_sp_dest_1, 
    dm.X_sp_orig_2, 
    dm.X_sp_dest_2, 
    dm.Z_sp_orig_2, 
    dm.Z_sp_dest_2, 
    Insample,
    dm.ID_1, 
    dm.ID_2,
    dm.ID_sp_dest_1, 
    dm.ID_sp_dest_2,
    X_orig,
    X_dest,
    Z_orig,
    s,
    s_sp_orig,
    s_sp_dest,
    s_m_cpl,
    s_f_cpl,
    insample,
    td,
    t2,
    t,
    sex,
    married,
    married_tminus1,
    x_orig_fhat,
    x_dest_fhat,
    xstate_orig_fhat,
    xstate_dest_fhat,
    y_orig_fhat,
    z_orig_fhat,
    zstate_orig_fhat,
    x_sp_orig_fhat,
    x_sp_dest_fhat,
    xstate_sp_orig_fhat,
    xstate_sp_dest_fhat,
    y_sp_orig_fhat,
    z_sp_orig_fhat,
    zstate_sp_orig_fhat,
    x_m_cpl,
    x_f_cpl,
    xstate_m_cpl,
    xstate_f_cpl,
    z_m_cpl,
    z_f_cpl,
    zstate_m_cpl,
    zstate_f_cpl,
    x,
    xstate,
    z,
    zstate,
    y,
    y_m_cpl,
    y_f_cpl,
    x_sp_orig,
    xstate_sp_orig,
    z_sp_orig,
    zstate_sp_orig,
    x_sp_dest,
    xstate_sp_dest,
    z_sp_dest,
    zstate_sp_dest,
    id,
    firstsampled,
    I_st1,
    I_st2_single,
    I_st2_cpl)

    return  data
end

struct Data
    Sex_1:: Union{Array,Matrix,Vector}
    X_1:: Union{Array,Matrix,Vector}
    Z_1:: Union{Array,Matrix,Vector}
    State_1:: Union{Array,Matrix,Vector}
    Zstate_1:: Union{Array,Matrix,Vector}
    Xstate_1:: Union{Array,Matrix,Vector}
    Y_1:: Union{Array,Matrix,Vector}
    Married_1:: Union{Array,Matrix,Vector}
    Sampled_1:: Union{Array,Matrix,Vector}
    Sex_2:: Union{Array,Matrix,Vector}
    X_2:: Union{Array,Matrix,Vector} 
    Z_2:: Union{Array,Matrix,Vector} 
    State_2:: Union{Array,Matrix,Vector} 
    Zstate_2:: Union{Array,Matrix,Vector} 
    Xstate_2:: Union{Array,Matrix,Vector} 
    Y_2:: Union{Array,Matrix,Vector}
    Married_2:: Union{Array,Matrix,Vector}
    Sampled_2:: Union{Array,Matrix,Vector}
    Zstate_sp_orig_1:: Union{Array,Matrix,Vector} 
    Zstate_sp_dest_1:: Union{Array,Matrix,Vector} 
    Xstate_sp_orig_1:: Union{Array,Matrix,Vector} 
    Xstate_sp_dest_1:: Union{Array,Matrix,Vector} 
    Y_sp_dest_1:: Union{Array,Matrix,Vector}
    Zstate_sp_orig_2:: Union{Array,Matrix,Vector} 
    Zstate_sp_dest_2:: Union{Array,Matrix,Vector} 
    Xstate_sp_orig_2:: Union{Array,Matrix,Vector} 
    Xstate_sp_dest_2:: Union{Array,Matrix,Vector} 
    Y_sp_dest_2:: Union{Array,Matrix,Vector}
    X_sp_orig_1:: Union{Array,Matrix,Vector} 
    X_sp_dest_1:: Union{Array,Matrix,Vector} 
    Z_sp_orig_1:: Union{Array,Matrix,Vector} 
    Z_sp_dest_1:: Union{Array,Matrix,Vector} 
    X_sp_orig_2:: Union{Array,Matrix,Vector} 
    X_sp_dest_2:: Union{Array,Matrix,Vector} 
    Z_sp_orig_2:: Union{Array,Matrix,Vector} 
    Z_sp_dest_2:: Union{Array,Matrix,Vector} 
    Insample:: Union{Array,Matrix,Vector}
    ID_1:: Union{Array,Matrix,Vector} 
    ID_2:: Union{Array,Matrix,Vector}
    ID_sp_dest_1:: Union{Array,Matrix,Vector} 
    ID_sp_dest_2:: Union{Array,Matrix,Vector}
    X_orig:: Union{Array,Matrix,Vector}
    X_dest:: Union{Array,Matrix,Vector}
    Z_orig:: Union{Array,Matrix,Vector}
    s:: Union{Array,Matrix,Vector} 
    s_sp_orig:: Union{Array,Matrix,Vector} 
    s_sp_dest:: Union{Array,Matrix,Vector} 
    s_m_cpl:: Union{Array,Matrix,Vector} 
    s_f_cpl:: Union{Array,Matrix,Vector} 
    insample:: Union{Array,Matrix,Vector} 
    td:: Union{Array,Matrix,Vector} 
    t2:: Union{Array,Matrix,Vector}
    t:: Union{Array,Matrix,Vector}
    sex:: Union{Array,Matrix,Vector}
    married:: Union{Array,Matrix,Vector}
    married_tminus1:: Union{Array,Matrix,Vector}
    x_orig_fhat:: Union{Array,Matrix,Vector}
    x_dest_fhat:: Union{Array,Matrix,Vector}
    xstate_orig_fhat:: Union{Array,Matrix,Vector}
    xstate_dest_fhat:: Union{Array,Matrix,Vector}
    y_orig_fhat:: Union{Array,Matrix,Vector}
    z_orig_fhat:: Union{Array,Matrix,Vector}
    zstate_orig_fhat:: Union{Array,Matrix,Vector}
    x_sp_orig_fhat:: Union{Array,Matrix,Vector}
    x_sp_dest_fhat:: Union{Array,Matrix,Vector}
    xstate_sp_orig_fhat:: Union{Array,Matrix,Vector}
    xstate_sp_dest_fhat:: Union{Array,Matrix,Vector}
    y_sp_orig_fhat:: Union{Array,Matrix,Vector}
    z_sp_orig_fhat:: Union{Array,Matrix,Vector}
    zstate_sp_orig_fhat:: Union{Array,Matrix,Vector}
    x_m_cpl:: Union{Array,Matrix,Vector}
    x_f_cpl:: Union{Array,Matrix,Vector}
    xstate_m_cpl:: Union{Array,Matrix,Vector}
    xstate_f_cpl:: Union{Array,Matrix,Vector}
    z_m_cpl:: Union{Array,Matrix,Vector}
    z_f_cpl:: Union{Array,Matrix,Vector}
    zstate_m_cpl:: Union{Array,Matrix,Vector}
    zstate_f_cpl:: Union{Array,Matrix,Vector}
    x:: Union{Array,Matrix,Vector}
    xstate:: Union{Array,Matrix,Vector}
    z:: Union{Array,Matrix,Vector}
    zstate:: Union{Array,Matrix,Vector}
    y:: Union{Array,Matrix,Vector}
    y_m_cpl:: Union{Array,Matrix,Vector}
    y_f_cpl:: Union{Array,Matrix,Vector}
    x_sp_orig:: Union{Array,Matrix,Vector}
    xstate_sp_orig:: Union{Array,Matrix,Vector}
    z_sp_orig:: Union{Array,Matrix,Vector}
    zstate_sp_orig:: Union{Array,Matrix,Vector}
    x_sp_dest:: Union{Array,Matrix,Vector}
    xstate_sp_dest:: Union{Array,Matrix,Vector}
    z_sp_dest:: Union{Array,Matrix,Vector}
    zstate_sp_dest:: Union{Array,Matrix,Vector}
    id:: Union{Array,Matrix,Vector}
    firstsampled:: Union{Array,Matrix,Vector}
    I_st1:: Union{Array,Matrix,Vector}
    I_st2_single:: Union{Array,Matrix,Vector}
    I_st2_cpl:: Union{Array,Matrix,Vector}
end


