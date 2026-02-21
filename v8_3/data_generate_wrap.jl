function f_data_generate_wrap(CCP_true_sim=[],load_=false)
    
    if load_==false
        if !isdefined(Main, :CCP_true_sim)
            CCP_true_sim=f_data_genbus4_solve()
        end
        dm = f_data_genbus4_simulate(CCP_true_sim,N_)
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
    # Spouses identifiers
    ###############################################
    # Notice that we always know what individual each observation pertains to based on their position in the vector
    ID_sp_dest=vcat(dm.ID_sp_dest_1,dm.ID_sp_dest_2)
    id_sp_dest=ID_sp_dest[:]

    ###########################################################
    # Indicator for when a person is in the sample or not
    ###########################################################
    # Each originally sampled agent faces an exogenous sample drop-out probability
    # Furthermore, spouses of originally sample agents drop out of the sample if they are not married to them anymore.

    Insample=vcat(dm.Sampled_1,dm.Sampled_2)
    insample=Insample[:]
    insample=repeat(insample,NTypes^2)
    insample!_Bool=insample.==0

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
    firstsampled_Bool=firstsampled.==1

    ###########################################################
    # Indicator for when a person is married or not
    ###########################################################
    Married=vcat(dm.Married_1,dm.Married_2)
    # married=repeat(coalesce.(Married[:], NaN), NTypes^2)
    married=repeat(Married[:],NTypes^2)

    Married_tminus1_1=hcat(Array{Union{Missing, Int64}}(missing, N),dm.Married_1[:,2:end])
    Married_tminus1_2=hcat(Array{Union{Missing, Int64}}(missing, N),dm.Married_2[:,2:end])
    Married_tminus1=vcat(Married_tminus1_1,Married_tminus1_2)
    # married_tminus1=repeat(coalesce.(Married_tminus1[:], NaN), NTypes^2)
    married_tminus1=repeat(Married_tminus1[:],NTypes^2)
    
    #####################################################################
    # Indicator for whether a household is composed of a single person
    # or a married couple when first sampled
    #####################################################################
    marriedhh=coalesce.(married,false).*firstsampled
    marriedhh=maximum(reshape(marriedhh,(2*N,T2,NTypes^2)),dims=2)
    MarriedHH=marriedhh[:,1,1]
    marriedhh=repeat(MarriedHH,1,T2,NTypes^2)[:]

    #####################################################################
    # Indicator for each period spent with original spouse from the 
    # time of first sample
    #####################################################################
    OriginalSpouse=vcat(dm.ID_sp_dest_1,dm.ID_sp_dest_2).==vcat(dm.ID_2,dm.ID_1)
    OriginalSpouse=repeat(OriginalSpouse,1,1,NTypes,NTypes_sp)
    OriginalSpouse=reshape(OriginalSpouse,(N,2,T2,NTypes,NTypes_sp)) # This is used in update_ptype
    OriginalSpouse=copy(OriginalSpouse)
    originalspouse=coalesce.(OriginalSpouse[:], false)

    ###########################################################
    # Sex
    ###########################################################
    Sex=vcat(dm.Sex_1,dm.Sex_2)
    sex=repeat(Sex,T2)[:]
    sex=repeat(sex,NTypes^2)

    ###########################################################
    # Spouse id
    ###########################################################
    id=vcat(dm.ID_sp_dest_1,dm.ID_sp_dest_2)[:]

    ###########################################################
    # Regressors for transitions estimation
    ###########################################################
    
    # Agents' transitions
    X_dest=vcat(dm.Phi_1[:,2:end,1],dm.Phi_2[:,2:end,1])
    X_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_dest)
    X_orig=vcat(dm.Phi_1[:,1:end-1,1],dm.Phi_2[:,1:end-1,1])
    X_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_orig)
    Phistate_dest=vcat(dm.Phistate_1[:,2:end],dm.Phistate_2[:,2:end])
    Phistate_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Phistate_dest)
    Phistate_orig=vcat(dm.Phistate_1[:,1:end-1],dm.Phistate_2[:,1:end-1])
    Phistate_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Phistate_orig)
    
    H_orig=vcat(repeat(dm.H_1,1,T2)[:,2:end],repeat(dm.H_2,1,T2)[:,2:end])
    H_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),H_orig)
    Hstate_orig=vcat(repeat(dm.Hstate_1,1,T2)[:,2:end],repeat(dm.Hstate_2,1,T2)[:,2:end])
    Hstate_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Hstate_orig)

    Y_orig=vcat(dm.Y_1[:,1:end-1],dm.Y_2[:,1:end-1])
    Y_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Y_orig)
    
    x_orig_fhat=repeat(X_orig[:],NTypes^2)
    x_dest_fhat=repeat(X_dest[:],NTypes^2)
    xstate_orig_fhat=repeat(Phistate_orig[:],NTypes^2)
    xstate_dest_fhat=repeat(Phistate_dest[:],NTypes^2)
    h_orig_fhat=repeat(H_orig[:],NTypes^2)
    hstate_orig_fhat=repeat(Hstate_orig[:],NTypes^2)
    y_orig_fhat=repeat(Y_orig[:],NTypes^2)

    # Remember that by construction, replacing the engine (y=1) instantaneously sets mileage to zero, so
    x_orig_fhat[y_orig_fhat.===1].=0
    xstate_orig_fhat[y_orig_fhat.===1].=1

    # Spouses' transitions
    # The matrix Phi_sp_dest is built using Phi_sp_orig_1 because in Phi_sp_orig_1 "orig" refers to the origin spouse, i.e. the spouse chosen at t-1, but when thinking of state transitions, Phi_sp_orig_1 contains the destination *state* of the spouse as they transition from t-1 to t.
    X_sp_dest=vcat(dm.Phi_sp_orig_1[:,2:end,1],dm.Phi_sp_orig_2[:,2:end,1])
    X_sp_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_sp_dest)
    X_sp_orig=vcat(dm.Phi_sp_dest_1[:,1:end-1,1],dm.Phi_sp_dest_2[:,1:end-1,1])
    X_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), X_sp_orig)
    Phistate_sp_dest=vcat(dm.Phistate_sp_orig_1[:,2:end],dm.Phistate_sp_orig_2[:,2:end])
    Phistate_sp_dest=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Phistate_sp_dest)
    Phistate_sp_orig=vcat(dm.Phistate_sp_dest_1[:,1:end-1],dm.Phistate_sp_dest_2[:,1:end-1])
    Phistate_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Phistate_sp_orig)
    
    H_sp_orig=vcat(dm.H_sp_orig_1[:,2:end],dm.H_sp_orig_2[:,2:end])
    H_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), H_sp_orig)
    Hstate_sp_orig=vcat(dm.Hstate_sp_orig_1[:,2:end],dm.Hstate_sp_orig_2[:,2:end])
    Hstate_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N), Hstate_sp_orig)

    Y_sp_orig=vcat(dm.Y_sp_dest_1[:,1:end-1],dm.Y_sp_dest_2[:,1:end-1])
    Y_sp_orig=hcat(Array{Union{Missing, Int64}}(missing, 2*N),Y_sp_orig)
    
    x_sp_orig_fhat=repeat(X_sp_orig[:],NTypes^2)
    x_sp_dest_fhat=repeat(X_sp_dest[:],NTypes^2)
    xstate_sp_orig_fhat=repeat(Phistate_sp_orig[:],NTypes^2)
    xstate_sp_dest_fhat=repeat(Phistate_sp_dest[:],NTypes^2)
    h_sp_orig_fhat=repeat(H_sp_orig[:],NTypes^2)
    hstate_sp_orig_fhat=repeat(Hstate_sp_orig[:],NTypes^2)
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
    x=vcat(dm.Phi_1[:,:,1],dm.Phi_2[:,:,1])[:] 
    temp_1=vcat(dm.Phistate_1,dm.Phistate_2)[:]
    temp_2=coalesce.(temp_1,1)
    xstate=Vector{Union{Missing,Int64}}(undef,size(x))
    xstate=ifelse.(temp_1.===missing,missing,phistate_list[temp_2,1])
    H_1_temp=ifelse.(dm.Sampled_1.==0,missing,repeat(dm.H_1,1,T2))
    H_2_temp=ifelse.(dm.Sampled_2.==0,missing,repeat(dm.H_2,1,T2))
    h=vcat(H_1_temp,H_2_temp)[:]
    Hstate_1_temp=ifelse.(dm.Sampled_1.==0,missing,repeat(dm.Hstate_1,1,T2))
    Hstate_2_temp=ifelse.(dm.Sampled_2.==0,missing,repeat(dm.Hstate_2,1,T2))
    hstate=vcat(Hstate_1_temp,Hstate_2_temp)[:]
    y=vcat(dm.Y_1,dm.Y_2)[:]

    x=repeat(x,NTypes^2)
    xstate=repeat(xstate,NTypes^2)
    h=repeat(h,NTypes^2)
    hstate=repeat(hstate,NTypes^2)
    y=repeat(y,NTypes^2)


    # Couple characteristics for choice probabilities as married
    ## X
    x_sp_temp=vcat(dm.Phi_sp_dest_1[:,:,1],dm.Phi_sp_dest_2[:,:,1])[:]
    x_sp_temp=repeat(x_sp_temp,NTypes^2)

    temp_1=vcat(dm.Phistate_sp_dest_1,dm.Phistate_sp_dest_2)[:]
    temp_2=coalesce.(temp_1,1)
    xstate_sp_temp=Vector{Union{Missing,Int64}}(undef,size(x_sp_temp))
    xstate_sp_temp=ifelse.(temp_1.===missing,missing,phistate_list[temp_2,1])
    
    xstate_sp_temp=repeat(xstate_sp_temp,NTypes^2)

    x_m_cpl= x .* (sex.==1) .+ x_sp_temp .* (sex.==2)
    x_f_cpl= x .* (sex.==2) .+ x_sp_temp .* (sex.==1)

    xstate_m_cpl= xstate .* (sex.==1) + xstate_sp_temp .* (sex.==2)
    xstate_f_cpl= xstate .* (sex.==2) + xstate_sp_temp .* (sex.==1)

    ## Z
    h_sp_temp=vcat(dm.H_sp_dest_1,dm.H_sp_dest_2)[:]
    h_sp_temp=repeat(h_sp_temp,NTypes^2)
    hstate_sp_temp=vcat(dm.Hstate_sp_dest_1,dm.Hstate_sp_dest_2)[:]
    hstate_sp_temp=repeat(hstate_sp_temp,NTypes^2)

    h_m_cpl= h .* (sex.==1) .+ h_sp_temp .* (sex.==2)
    h_f_cpl= h .* (sex.==2) .+ h_sp_temp .* (sex.==1)

    hstate_m_cpl= hstate .* (sex.==1) + hstate_sp_temp .* (sex.==2)
    hstate_f_cpl= hstate .* (sex.==2) + hstate_sp_temp .* (sex.==1)

    ## Y
    y_sp_temp=vcat(dm.Y_sp_dest_1,dm.Y_sp_dest_2)[:]
    y_sp_temp=repeat(y_sp_temp,NTypes^2)

    y_m_cpl= y .* (sex.==1) .+ y_sp_temp .* (sex.==2)
    y_f_cpl= y .* (sex.==2) .+ y_sp_temp .* (sex.==1)

    ###########################################################
    # Regressors for first-stage choice probabilities
    ###########################################################
    phi=Matrix{Union{Missing,Float64}}(undef,2*N*T2,phidim)
    for dim in 1:phidim
        phi[:,dim]=vcat(dm.Phi_1[:,:,dim],dm.Phi_2[:,:,dim])[:]
    end
    phi=repeat(phi,NTypes^2)
    phistate=vcat(dm.Phistate_1,dm.Phistate_2)[:]
    phistate=repeat(phistate,NTypes^2)
    
    phi_sp_orig=Matrix{Union{Missing,Float64}}(undef,2*N*T2,phidim)
    for dim in 1:phidim
        phi_sp_orig[:,dim]=vcat(dm.Phi_sp_orig_1[:,:,dim],dm.Phi_sp_orig_2[:,:,dim])[:]
    end
    phi_sp_orig=repeat(phi_sp_orig,NTypes^2)
    phistate_sp_orig=vcat(dm.Phistate_sp_orig_1,dm.Phistate_sp_orig_2)[:]
    phistate_sp_orig=repeat(phistate_sp_orig,NTypes^2)
    if COND_ORIG_SPOUSE==0
        phi_sp_orig=zeros(Union{Missing,Float64},size(phi_sp_orig))
        phistate_sp_orig=ones(Union{Missing,Int64},size(phistate_sp_orig))
    end

    phi_sp_dest=Matrix{Union{Missing,Float64}}(undef,2*N*T2,phidim)
    for dim in 1:phidim
        phi_sp_dest[:,dim]=vcat(dm.Phi_sp_dest_1[:,:,dim],dm.Phi_sp_dest_2[:,:,dim])[:]
    end
    phi_sp_dest=repeat(phi_sp_dest,NTypes^2)
    phistate_sp_dest=vcat(dm.Phistate_sp_dest_1,dm.Phistate_sp_dest_2)[:]
    phistate_sp_dest=repeat(phistate_sp_dest,NTypes^2)
    

    h_sp_orig=vcat(dm.H_sp_orig_1,dm.H_sp_orig_2)[:]
    h_sp_orig=repeat(h_sp_orig,NTypes^2)
    hstate_sp_orig=vcat(dm.Hstate_sp_orig_1,dm.Hstate_sp_orig_2)[:]
    hstate_sp_orig=repeat(hstate_sp_orig,NTypes^2)
    if COND_ORIG_SPOUSE==0
        h_sp_orig=zeros(Union{Missing,Float64},size(h_sp_orig))
        hstate_sp_orig=ones(Union{Missing,Int64},size(hstate_sp_orig))
    end
    h_sp_dest=vcat(dm.H_sp_dest_1,dm.H_sp_dest_2)[:]
    h_sp_dest=repeat(h_sp_dest,NTypes^2)
    hstate_sp_dest=vcat(dm.Hstate_sp_dest_1,dm.Hstate_sp_dest_2)[:]
    hstate_sp_dest=repeat(hstate_sp_dest,NTypes^2)

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
    s=ifelse.((insample.==0),missing,s)

    ############################
    # Used in first step CCPs
    ############################
    s_sp_orig=ones(Int, size(s))
    s_sp_dest=repeat(s_t_sp_dest , T2,1)[:]
    s_sp_dest=ifelse.(((married.===0) .| (insample.==0)), missing, s_sp_dest)
    
    ############################
    # Used in second step CCPs
    ############################
    s_m_cpl = s.*(sex.===1) .+ s_sp_dest.*(sex.===2)
    # s_m_cpl = ifelse.(((married.===0) .| (insample.==0))[:], NaN, s_m_cpl) 
    s_m_cpl=ifelse.(((married.===0) .| (insample.==0))[:], missing, s_m_cpl)
    s_f_cpl = s.*(sex.===2) .+ s_sp_dest.*(sex.===1)
    # s_f_cpl = ifelse.(((married.===0) .| (insample.==0))[:], NaN, s_f_cpl) 
    s_f_cpl = ifelse.(((married.===0) .| (insample.==0))[:], missing, s_f_cpl)

    ############################################### 
    # True value of unobserved state 
    ############################################### 
    State=vcat(dm.State_1,dm.State_2)
    State=repeat(State,1, T2)
    s_true=repeat(State[:], NTypes^2)

    State_sp_dest=vcat(dm.State_sp_dest_1,dm.State_sp_dest_2)
    s_sp_dest_true=repeat(State_sp_dest[:], NTypes^2)

    s_m_cpl_true=s_true.*(sex.===1) .+ s_sp_dest_true.*(sex.===2)
    s_m_cpl_true=ifelse.(((married.===0) .| (insample.==0)),missing,s_m_cpl_true)
    s_f_cpl_true=s_true.*(sex.===2) .+ s_sp_dest_true.*(sex.===1)
    s_f_cpl_true=ifelse.(((married.===0) .| (insample.==0)),missing,s_f_cpl_true)
    
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
    I_st1=Array{Union{Missing, Int}}(missing, size(x)[1],phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest)
    for obs in 1:size(x)[1]
        for phistate_sp in 1:phibin_sp_dest
            for hstate_sp in 1:hbin_sp_dest
                for s_sp in 1:NTypes_sp_dest
                    # Triple inequality because these vectors might contain missing values
                    if ismissing(phistate_sp_dest[obs].*hstate_sp_dest[obs].*s_sp_dest[obs])
                        continue
                    end
                    I_st1[obs,phistate_sp,hstate_sp,s_sp]=(phistate_sp_dest[obs]===phistate_sp_dest)*(hstate_sp_dest[obs]===hstate_sp_dest)*(s_sp_dest[obs]===s_sp_dest)
                end
            end
        end
    end
    # I_st1 will need to be comformable to element-wise subtraction with the CCP arrays CCP_struct.p_st2_f and CCP_struct.p_st2_m, so I need to add some placeholder 

    # Second-stage: one set of indicators for each choice the agent can choose as single, and one set of indicators for each choice a couple can choose
    I_st2_single=Array{Union{Missing, Int}}(missing, size(x)[1],J+1)
    for obs in 1:size(x)[1]
        for a in 1:J+1
            # Triple inequality because these vectors might contain missing values
            if ismissing(y[obs])
                continue
            end
            I_st2_single[obs,a]=(y[obs]===a)
        end
    end

    I_st2_cpl=Array{Union{Missing, Int}}(missing, size(x)[1],J_m+1,J_f+1)
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


    ########################################################################################################
    #   Coalesced versions of regression and indexes
    ########################################################################################################
    phistate_sp_dest_=coalesce.(phistate_sp_dest,1)
    hstate_sp_dest_=coalesce.(hstate_sp_dest,1)
    s_sp_dest_=coalesce.(s_sp_dest,1)
    phistate_=coalesce.(phistate,1)
    xstate_=coalesce.(xstate,1)
    hstate_=coalesce.(hstate,1)
    s_=coalesce.(s,1)
    phistate_sp_orig_=coalesce.(phistate_sp_orig,1)
    hstate_sp_orig_=coalesce.(hstate_sp_orig,1)
    s_sp_orig_=coalesce.(s_sp_orig,1)
    xstate_m_cpl_=coalesce.(xstate_m_cpl,1)
    hstate_m_cpl_=coalesce.(hstate_m_cpl,1)
    s_m_cpl_=coalesce.(s_m_cpl,1)
    xstate_f_cpl_=coalesce.(xstate_f_cpl,1)
    hstate_f_cpl_=coalesce.(hstate_f_cpl,1)
    s_f_cpl_=coalesce.(s_f_cpl,1)
    y_=coalesce.(y,1)
    y_m_cpl_=coalesce.(y_m_cpl,1)
    y_f_cpl_=coalesce.(y_f_cpl,1)
    married_=coalesce.(married,0)
    sex_=coalesce.(sex,0)
    married_tminus1_=coalesce.(married_tminus1,0)
    I_st1_=coalesce.(I_st1,0)
    I_st2_single_=coalesce.(I_st2_single,0)
    I_st2_cpl_=coalesce.(I_st2_cpl,0)
    Sex_=coalesce.(Sex,0)
    OriginalSpouse_=coalesce.(OriginalSpouse,false)

    ########################################################################################################
    #   Indicators, Cartesian indexes, and Boolean arrays used in the update of posterior type probabilities
    ########################################################################################################
    same_sp = comp.(eachcol(ID_sp_dest), Ref(ID_sp_dest))
    same_sp=cat(same_sp..., dims=3)

    First_t_with_sp= [[findfirst(Ref(ID_sp_dest[i,j]) .=== ID_sp_dest[i, :]).==j for j in 1:size(ID_sp_dest,2)] for i in 1:size(ID_sp_dest,1)]
    First_t_with_sp=vcat(transpose.(First_t_with_sp)...)
    First_t_with_sp=repeat(First_t_with_sp,1,1,NTypes,NTypes_sp)
    First_t_with_sp=reshape(First_t_with_sp,(N,2,T2,NTypes,NTypes_sp))

    NotMarried= coalesce.(1 .-Married,0)
    NotMarried=reshape(NotMarried,(N,2,T2))

    NotMarried_m=NotMarried[:,[1],:].*(Sex_[1:N,1].==1) .+ NotMarried[:,[2],:].*(Sex_[1:N,1].==2)
    NotMarried_f=NotMarried[:,[1],:].*(Sex_[1:N,1].==2) .+ NotMarried[:,[2],:].*(Sex_[1:N,1].==1)
    NotMarried_sort=cat(NotMarried_m,NotMarried_f,dims=2)
    NotMarried_sort=repeat(NotMarried_sort,1,1,1,NTypes_m,NTypes_f)    

    OriginalSpouse_m=OriginalSpouse_[:,[1],:,:,:].*(Sex_[1:N,1].==1) .+ OriginalSpouse_[:,[2],:,:,:].*(Sex_[1:N,1].==2)
    OriginalSpouse_f=OriginalSpouse_[:,[1],:,:,:].*(Sex_[1:N,1].==2) .+ OriginalSpouse_[:,[2],:,:,:].*(Sex_[1:N,1].==1)
    OriginalSpouse_sort=cat(OriginalSpouse_m,OriginalSpouse_f, dims=2)

    insample_temp=reshape(insample,N,2,T2,NTypes,NTypes_sp)
    insample_sort_m=insample_temp[:,[1],:,:,:].*(Sex_[1:N,1].==1) .+ insample_temp[:,[2],:,:,:].*(Sex_[1:N,1].==2)
    insample_sort_f=insample_temp[:,[1],:,:,:].*(Sex_[1:N,1].==2) .+ insample_temp[:,[2],:,:,:].*(Sex_[1:N,1].==1)
    insample_sort=cat(insample_sort_m,insample_sort_f, dims=2)

    idx_NotMarried_sort=LinearIndices((N,2,T2,NTypes,NTypes))[findall(x -> x == 1, NotMarried_sort)]
    idx_OriginalSpouse_sort=LinearIndices((N,2,T2,NTypes,NTypes))[findall(x -> x == 1, OriginalSpouse_sort)]
    idx_insample_sort=LinearIndices((N,2,T2,NTypes,NTypes))[findall(x -> x == 0, insample_sort)]

    ################################################################
    # Cartesian indexes indicating position in state space
    ################################################################

    # Indexes used to calculate the CCP corresponding to each observed choice. Used in likeCCP
    idx_p_st1_m_singleb4=CartesianIndex.(phistate_sp_dest_,hstate_sp_dest_,s_sp_dest_,phistate_,hstate_,s_,t)
    idx_p_st1_f_singleb4=CartesianIndex.(phistate_sp_dest_,hstate_sp_dest_,s_sp_dest_,phistate_,hstate_,s_,t)

    idx_p0_st1_m_singleb4=CartesianIndex.(phistate_,hstate_,s_,t)
    idx_p0_st1_f_singleb4=CartesianIndex.(phistate_,hstate_,s_,t)

    idx_p_st1_m_marriedb4=CartesianIndex.(phistate_sp_dest_,hstate_sp_dest_,s_sp_dest_,phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t)
    idx_p_st1_f_marriedb4=CartesianIndex.(phistate_sp_dest_,hstate_sp_dest_,s_sp_dest_,phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t)

    idx_p0_st1_m_marriedb4=CartesianIndex.(phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t)   
    idx_p0_st1_f_marriedb4=CartesianIndex.(phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t)   

    idx_p_st2_m_single=CartesianIndex.(xstate_,hstate_,s_,y_,t)
    idx_p_st2_f_single=CartesianIndex.(xstate_,hstate_,s_,y_,t)

    idx_p_st2_cpl=CartesianIndex.(xstate_m_cpl_,hstate_m_cpl_,s_m_cpl_,xstate_f_cpl_,hstate_f_cpl_,s_f_cpl_,y_m_cpl_,y_f_cpl_,t)

    # Reduced version of the above indexes, they include only the observations relevant to each choice probability.
    # This speeds up likeCCP
    idx_p_st1_m_singleb4_reduced_obs=(sex.===1).&&(firstsampled.==0).&&(married_tminus1.===0).&&(married.===1).&&(insample.==1)
    if sum(idx_p_st1_m_singleb4_reduced_obs)>1
        idx_p_st1_m_singleb4_reduced=@view idx_p_st1_m_singleb4[idx_p_st1_m_singleb4_reduced_obs]
    else
        idx_p_st1_m_singleb4_reduced=zeros(0)
    end
    idx_p_st1_f_singleb4_reduced_obs=(sex.===2).&&(firstsampled.==0).&&(married_tminus1.===0).&&(married.===1).&&(insample.==1)
    if sum(idx_p_st1_f_singleb4_reduced_obs)>1
        idx_p_st1_f_singleb4_reduced=@view idx_p_st1_f_singleb4[idx_p_st1_f_singleb4_reduced_obs]
    else
        idx_p_st1_f_singleb4_reduced=zeros(0)
    end
    idx_p0_st1_m_singleb4_reduced_obs=(sex.===1).&&(firstsampled.==0).&&(married_tminus1.===0).&&(married.===0).&&(insample.==1)
    if sum(idx_p0_st1_m_singleb4_reduced_obs)>1
        idx_p0_st1_m_singleb4_reduced=@view idx_p0_st1_m_singleb4[idx_p0_st1_m_singleb4_reduced_obs]
    else
        idx_p0_st1_m_singleb4_reduced=zeros(0)
    end
    idx_p0_st1_f_singleb4_reduced_obs=(sex.===2).&&(firstsampled.==0).&&(married_tminus1.===0).&&(married.===0).&&(insample.==1)
    if sum(idx_p0_st1_f_singleb4_reduced_obs)>1
        idx_p0_st1_f_singleb4_reduced=@view idx_p0_st1_f_singleb4[idx_p0_st1_f_singleb4_reduced_obs]
    else
        idx_p0_st1_f_singleb4_reduced=zeros(0)
    end
    idx_p_st1_m_marriedb4_reduced_obs=(sex.===1).&&(firstsampled.==0).&&(married_tminus1.===1).&&(married.===1).&&(insample.==1)
    if sum(idx_p_st1_m_marriedb4_reduced_obs)>1
        idx_p_st1_m_marriedb4_reduced=@view idx_p_st1_m_marriedb4[idx_p_st1_m_marriedb4_reduced_obs]
    else
        idx_p_st1_m_marriedb4_reduced=zeros(0)
    end
    idx_p_st1_f_marriedb4_reduced_obs=(sex.===2).&&(firstsampled.==0).&&(married_tminus1.===1).&&(married.===1).&&(insample.==1)
    if sum(idx_p_st1_f_marriedb4_reduced_obs)>1
        idx_p_st1_f_marriedb4_reduced=@view idx_p_st1_f_marriedb4[idx_p_st1_f_marriedb4_reduced_obs]
    else
        idx_p_st1_f_marriedb4_reduced=zeros(0)
    end
    idx_p0_st1_m_marriedb4_reduced_obs=(sex.===1).&&(firstsampled.==0).&&(married_tminus1.===1).&&(married.===0).&&(insample.==1)
    if sum(idx_p0_st1_m_marriedb4_reduced_obs)>1
        idx_p0_st1_m_marriedb4_reduced=@view idx_p0_st1_m_marriedb4[idx_p0_st1_m_marriedb4_reduced_obs]
    else
        idx_p0_st1_m_marriedb4_reduced=zeros(0)
    end
    idx_p0_st1_f_marriedb4_reduced_obs=(sex.===2).&&(firstsampled.==0).&&(married_tminus1.===1).&&(married.===0).&&(insample.==1)
    if sum(idx_p0_st1_f_marriedb4_reduced_obs)>1
        idx_p0_st1_f_marriedb4_reduced=@view idx_p0_st1_f_marriedb4[idx_p0_st1_f_marriedb4_reduced_obs]
    else
        idx_p0_st1_f_marriedb4_reduced=zeros(0)
    end
    idx_p_st2_m_single_reduced_obs=(sex.===1).&&(married.===0).&&(insample.==1)
    if sum(idx_p_st2_m_single_reduced_obs)>1
        idx_p_st2_m_single_reduced=@view idx_p_st2_m_single[idx_p_st2_m_single_reduced_obs]
    else
        idx_p_st2_m_single_reduced=zeros(0)
    end
    idx_p_st2_f_single_reduced_obs=(sex.===2).&&(married.===0).&&(insample.==1)
    if sum(idx_p_st2_f_single_reduced_obs)>1
        idx_p_st2_f_single_reduced=@view idx_p_st2_f_single[idx_p_st2_f_single_reduced_obs]
    else
        idx_p_st2_f_single_reduced=zeros(0)
    end
    idx_p_st2_cpl_reduced_obs=(married.===1).&&(insample.==1)
    if sum(idx_p_st2_cpl_reduced_obs)>1
        idx_p_st2_cpl_reduced=@view idx_p_st2_cpl[idx_p_st2_cpl_reduced_obs]
        originalspouse_reduced=@view originalspouse[idx_p_st2_cpl_reduced_obs]
    else
        idx_p_st2_cpl_reduced=zeros(0)
        originalspouse_reduced=zeros(0)
    end
    
    # Partial indexes used to calculate the CCP of each available choice at each observed state in the data. Used in resCCP
    partidx_p_st1_m_singleb4_data=CartesianIndex.(phistate_,hstate_,s_,t)
    partidx_p_st1_m_marriedb4_data=CartesianIndex.(phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t) 
    partidx_p_st1_f_singleb4_data=CartesianIndex.(phistate_,hstate_,s_,t)
    partidx_p_st1_f_marriedb4_data=CartesianIndex.(phistate_,hstate_,s_,phistate_sp_orig_,hstate_sp_orig_,s_sp_orig_,t) 
    partidx_p_st2_m_data=CartesianIndex.(xstate_,hstate_,s_,t)
    partidx_p_st2_f_data=CartesianIndex.(xstate_,hstate_,s_,t)
    partidx_p_st2_cpl_data=CartesianIndex.(xstate_m_cpl_,hstate_m_cpl_,s_m_cpl_,xstate_f_cpl_,hstate_f_cpl_,s_f_cpl_,t)
    
    # Indexes used to create the gradient of the CCPs with respect to reduced form parameters
    idx_p_st1_m_singleb4_choicesonly=  map(index -> (index[1], index[2], index[3]), idx_p_st1_m_singleb4)
    idx_p_st1_m_singleb4_choicesonly=CartesianIndex.(idx_p_st1_m_singleb4_choicesonly)
    idx_p_st1_m_marriedb4_choicesonly=  map(index -> (index[1], index[2], index[3]), idx_p_st1_m_marriedb4)
    idx_p_st1_m_marriedb4_choicesonly=CartesianIndex.(idx_p_st1_m_marriedb4_choicesonly)
    # Do not transform the following indexes to linear index. Used in reducedform_state_original
    idx_st1_m_singleb4_astuples = [Tuple(idx) for idx in  idx_p_st1_m_singleb4] 
    idx_st1_m_singleb4_astuples=CartesianIndex.(idx_st1_m_singleb4_astuples)
    idx_st1_m_marriedb4_astuples = [Tuple(idx) for idx in  idx_p_st1_m_marriedb4]
    idx_st1_m_marriedb4_astuples=CartesianIndex.(idx_st1_m_marriedb4_astuples)

    # Transform all indexes from cartesian into linear indexes. It's generally more memory efficient and also it works with CUDA
    idx_p_st1_m_singleb4=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))[idx_p_st1_m_singleb4]
    idx_p_st1_f_singleb4=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))[idx_p_st1_f_singleb4]
    idx_p0_st1_m_singleb4=LinearIndices((phibin_m,hbin_m,NTypes_m,T2))[idx_p0_st1_m_singleb4]
    idx_p0_st1_f_singleb4=LinearIndices((phibin_f,hbin_f,NTypes_f,T2))[idx_p0_st1_f_singleb4]
    idx_p_st1_m_marriedb4=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p_st1_m_marriedb4]
    idx_p_st1_f_marriedb4=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p_st1_f_marriedb4]
    idx_p0_st1_m_marriedb4=LinearIndices((phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p0_st1_m_marriedb4]  
    idx_p0_st1_f_marriedb4=LinearIndices((phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p0_st1_f_marriedb4]
    idx_p_st2_m_single=LinearIndices((xbin_m,hbin_m,NTypes_m,J_m+1,T2))[idx_p_st2_m_single]
    idx_p_st2_f_single=LinearIndices((xbin_f,hbin_f,NTypes_f,J_f+1,T2))[idx_p_st2_f_single]
    idx_p_st2_cpl=LinearIndices((xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))[idx_p_st2_cpl]
    
    idx_p_st1_m_singleb4_reduced=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))[idx_p_st1_m_singleb4_reduced]
    idx_p_st1_f_singleb4_reduced=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2))[idx_p_st1_f_singleb4_reduced]
    idx_p0_st1_m_singleb4_reduced=LinearIndices((phibin_m,hbin_m,NTypes_m,T2))[idx_p0_st1_m_singleb4_reduced]
    idx_p0_st1_f_singleb4_reduced=LinearIndices((phibin_f,hbin_f,NTypes_f,T2))[idx_p0_st1_f_singleb4_reduced]
    idx_p_st1_m_marriedb4_reduced=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p_st1_m_marriedb4_reduced]
    idx_p_st1_f_marriedb4_reduced=LinearIndices((phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p_st1_f_marriedb4_reduced]
    idx_p0_st1_m_marriedb4_reduced=LinearIndices((phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p0_st1_m_marriedb4_reduced]  
    idx_p0_st1_f_marriedb4_reduced=LinearIndices((phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[idx_p0_st1_f_marriedb4_reduced]
    idx_p_st2_m_single_reduced=LinearIndices((xbin_m,hbin_m,NTypes_m,J_m+1,T2))[idx_p_st2_m_single_reduced]
    idx_p_st2_f_single_reduced=LinearIndices((xbin_f,hbin_f,NTypes_f,J_f+1,T2))[idx_p_st2_f_single_reduced]
    idx_p_st2_cpl_reduced=LinearIndices((xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2))[idx_p_st2_cpl_reduced]

    # partidx_p_st1_m_singleb4_data=LinearIndices((phibin_m,hbin_m,NTypes_m,T2))[partidx_p_st1_m_singleb4_data]
    # partidx_p_st1_m_marriedb4_data=LinearIndices((phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[partidx_p_st1_m_marriedb4_data]
    # partidx_p_st1_f_singleb4_data=LinearIndices((phibin_f,hbin_f,NTypes_f,T2))[partidx_p_st1_f_singleb4_data]
    # partidx_p_st1_f_marriedb4_data=LinearIndices((phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))[partidx_p_st1_f_marriedb4_data]
    # partidx_p_st2_m_data=LinearIndices((xbin_m,hbin_m,NTypes_m,T2))[partidx_p_st2_m_data]
    # partidx_p_st2_f_data=LinearIndices((xbin_f,hbin_f,NTypes_f,T2))[partidx_p_st2_f_data]
    # partidx_p_st2_cpl_data=LinearIndices((xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))[partidx_p_st2_cpl_data]
    idx_p_st1_m_singleb4_choicesonly=LinearIndices((phibin_f,hbin_f,NTypes_f))[idx_p_st1_m_singleb4_choicesonly]
    idx_p_st1_m_marriedb4_choicesonly=LinearIndices((phibin_f,hbin_f,NTypes_f))[idx_p_st1_m_marriedb4_choicesonly]

    # Convert logical indexes into linear indexes
    idx_p_st1_m_singleb4_reduced_obs=findall(idx_p_st1_m_singleb4_reduced_obs.==1)
    idx_p_st1_f_singleb4_reduced_obs=findall(idx_p_st1_f_singleb4_reduced_obs.==1)
    idx_p0_st1_m_singleb4_reduced_obs=findall(idx_p0_st1_m_singleb4_reduced_obs.==1)
    idx_p0_st1_f_singleb4_reduced_obs=findall(idx_p0_st1_f_singleb4_reduced_obs.==1)
    idx_p_st1_m_marriedb4_reduced_obs=findall(idx_p_st1_m_marriedb4_reduced_obs.==1)
    idx_p_st1_f_marriedb4_reduced_obs=findall(idx_p_st1_f_marriedb4_reduced_obs.==1)
    idx_p0_st1_m_marriedb4_reduced_obs=findall(idx_p0_st1_m_marriedb4_reduced_obs.==1)
    idx_p0_st1_f_marriedb4_reduced_obs=findall(idx_p0_st1_f_marriedb4_reduced_obs.==1)
    idx_p_st2_m_single_reduced_obs=findall(idx_p_st2_m_single_reduced_obs.==1)
    idx_p_st2_f_single_reduced_obs=findall(idx_p_st2_f_single_reduced_obs.==1)
    idx_p_st2_cpl_reduced_obs=findall(idx_p_st2_cpl_reduced_obs.==1)

    #######################################################################################
    # Compress all vectors of indexes and binary, regardless of level of precision used.
    #######################################################################################

    same_sp=Bool.(same_sp)
    First_t_with_sp=Bool.(First_t_with_sp)
    NotMarried=Bool.(NotMarried)
    NotMarried_sort=Bool.(NotMarried_sort)
    OriginalSpouse_sort=Bool.(OriginalSpouse_sort)
    insample_sort=Bool.(insample_sort)
    idx_NotMarried_sort=compress_intarray(idx_NotMarried_sort)
    idx_OriginalSpouse_sort=compress_intarray(idx_OriginalSpouse_sort)
    idx_insample_sort=compress_intarray(idx_insample_sort)
    
    idx_p_st1_m_singleb4=compress_intarray(idx_p_st1_m_singleb4)
    idx_p_st1_f_singleb4=compress_intarray(idx_p_st1_f_singleb4)
    idx_p0_st1_m_singleb4=compress_intarray(idx_p0_st1_m_singleb4)
    idx_p0_st1_f_singleb4=compress_intarray(idx_p0_st1_f_singleb4)
    idx_p_st1_m_marriedb4=compress_intarray(idx_p_st1_m_marriedb4)
    idx_p_st1_f_marriedb4=compress_intarray(idx_p_st1_f_marriedb4)
    idx_p0_st1_m_marriedb4=compress_intarray(idx_p0_st1_m_marriedb4)
    idx_p0_st1_f_marriedb4=compress_intarray(idx_p0_st1_f_marriedb4)
    idx_p_st2_m_single=compress_intarray(idx_p_st2_m_single)
    idx_p_st2_f_single=compress_intarray(idx_p_st2_f_single)
    idx_p_st2_cpl=compress_intarray(idx_p_st2_cpl)
    
    idx_p_st1_m_singleb4_reduced=compress_intarray(idx_p_st1_m_singleb4_reduced)
    idx_p_st1_f_singleb4_reduced=compress_intarray(idx_p_st1_f_singleb4_reduced)
    idx_p0_st1_m_singleb4_reduced=compress_intarray(idx_p0_st1_m_singleb4_reduced)
    idx_p0_st1_f_singleb4_reduced=compress_intarray(idx_p0_st1_f_singleb4_reduced)
    idx_p_st1_m_marriedb4_reduced=compress_intarray(idx_p_st1_m_marriedb4_reduced)
    idx_p_st1_f_marriedb4_reduced=compress_intarray(idx_p_st1_f_marriedb4_reduced)
    idx_p0_st1_m_marriedb4_reduced=compress_intarray(idx_p0_st1_m_marriedb4_reduced)
    idx_p0_st1_f_marriedb4_reduced=compress_intarray(idx_p0_st1_f_marriedb4_reduced)
    idx_p_st2_m_single_reduced=compress_intarray(idx_p_st2_m_single_reduced)
    idx_p_st2_f_single_reduced=compress_intarray(idx_p_st2_f_single_reduced)
    idx_p_st2_cpl_reduced=compress_intarray(idx_p_st2_cpl_reduced)

    idx_p_st1_m_singleb4_choicesonly=compress_intarray(idx_p_st1_m_singleb4_choicesonly)
    idx_p_st1_m_marriedb4_choicesonly=compress_intarray(idx_p_st1_m_marriedb4_choicesonly)

    idx_p_st1_m_singleb4_reduced_obs=compress_intarray(idx_p_st1_m_singleb4_reduced_obs)
    idx_p_st1_f_singleb4_reduced_obs=compress_intarray(idx_p_st1_f_singleb4_reduced_obs)
    idx_p0_st1_m_singleb4_reduced_obs=compress_intarray(idx_p0_st1_m_singleb4_reduced_obs)
    idx_p0_st1_f_singleb4_reduced_obs=compress_intarray(idx_p0_st1_f_singleb4_reduced_obs)
    idx_p_st1_m_marriedb4_reduced_obs=compress_intarray(idx_p_st1_m_marriedb4_reduced_obs)
    idx_p_st1_f_marriedb4_reduced_obs=compress_intarray(idx_p_st1_f_marriedb4_reduced_obs)
    idx_p0_st1_m_marriedb4_reduced_obs=compress_intarray(idx_p0_st1_m_marriedb4_reduced_obs)
    idx_p0_st1_f_marriedb4_reduced_obs=compress_intarray(idx_p0_st1_f_marriedb4_reduced_obs)
    idx_p_st2_m_single_reduced_obs=compress_intarray(idx_p_st2_m_single_reduced_obs)
    idx_p_st2_f_single_reduced_obs=compress_intarray(idx_p_st2_f_single_reduced_obs)
    idx_p_st2_cpl_reduced_obs=compress_intarray(idx_p_st2_cpl_reduced_obs)

    MarriedHH=Bool.(MarriedHH)
    t=compress_intarray(t)
    phistate_sp_orig=compress_intarray(phistate_sp_orig)
    phistate_sp_orig_=compress_intarray(phistate_sp_orig_)
    hstate_sp_orig=compress_intarray(hstate_sp_orig)
    hstate_sp_orig_=compress_intarray(hstate_sp_orig_)
    s_sp_orig=compress_intarray(s_sp_orig)
    s_sp_orig_=compress_intarray(s_sp_orig_)
    s_sp_dest=compress_intarray(s_sp_dest)
    phistate=compress_intarray(phistate)
    phistate_=compress_intarray(phistate_)
    hstate=compress_intarray(hstate)
    hstate_=compress_intarray(hstate_)
    s=compress_intarray(s)
    s_=compress_intarray(s_)
    xstate=compress_intarray(xstate)
    xstate_=compress_intarray(xstate_)
    xstate_m_cpl=compress_intarray(xstate_m_cpl)
    xstate_m_cpl_=compress_intarray(xstate_m_cpl_)
    hstate_m_cpl=compress_intarray(hstate_m_cpl)
    hstate_m_cpl_=compress_intarray(hstate_m_cpl_)
    s_m_cpl=compress_intarray(s_m_cpl)
    s_m_cpl_=compress_intarray(s_m_cpl_)
    xstate_f_cpl=compress_intarray(xstate_f_cpl)
    xstate_f_cpl_=compress_intarray(xstate_f_cpl_)
    hstate_f_cpl=compress_intarray(hstate_f_cpl)
    hstate_f_cpl_=compress_intarray(hstate_f_cpl_)
    s_f_cpl=compress_intarray(s_f_cpl)
    s_f_cpl_=compress_intarray(s_f_cpl_)
    s_true=compress_intarray(s_true)
    s_m_cpl_true=compress_intarray(s_m_cpl_true)
    s_f_cpl_true=compress_intarray(s_f_cpl_true)
    insample=Bool.(insample)
    sex=compress_intarray(sex)
    s_true=compress_intarray(s_true)
    s=compress_intarray(s)
    s_sp_dest_true=compress_intarray(s_sp_dest_true)
    s_sp_dest=compress_intarray(s_sp_dest)

    xstate_orig_fhat=compress_intarray(xstate_orig_fhat)
    xstate_dest_fhat=compress_intarray(xstate_dest_fhat)
    hstate_orig_fhat=compress_intarray(hstate_orig_fhat)
    xstate_sp_orig_fhat=compress_intarray(xstate_sp_orig_fhat)
    xstate_sp_dest_fhat=compress_intarray(xstate_sp_dest_fhat)
    hstate_sp_orig_fhat=compress_intarray(hstate_sp_orig_fhat)
    firstsampled=Bool.(firstsampled)
    I_st1=map(x -> ismissing(x) ? missing : Bool(x), I_st1)
    I_st2_single=map(x -> ismissing(x) ? missing : Bool(x), I_st2_single)
    I_st2_cpl=map(x -> ismissing(x) ? missing : Bool(x), I_st2_cpl)
    I_st1_=Bool.(I_st1_)
    I_st2_single_=Bool.(I_st2_single_)
    I_st2_cpl_=Bool.(I_st2_cpl_)
    marriedhh=Bool.(marriedhh)
    married_tminus1_=Bool.(married_tminus1_)
    sex_=compress_intarray(sex_)
    married_=Bool.(married_)

    married = map(x -> ismissing(x) ? missing : Bool(x), married)
    married_tminus1= map(x -> ismissing(x) ? missing : Bool(x), married_tminus1)
    Married=map(x -> ismissing(x) ? missing : Bool(x), Married)

    
    
    data=Data(
        dm.Adj,
        dm.Sex_1,
        dm.Phi_1,
        dm.H_1,
        dm.State_1,
        dm.Hstate_1,
        dm.Phistate_1,
        dm.Y_1,
        dm.Married_1,
        dm.Sampled_1,
        dm.Sex_2,
        dm.Phi_2, 
        dm.H_2, 
        dm.State_2, 
        dm.Hstate_2, 
        dm.Phistate_2, 
        dm.Y_2,
        dm.Married_2,
        dm.Sampled_2,
        dm.Hstate_sp_orig_1, 
        dm.Hstate_sp_dest_1, 
        dm.Phistate_sp_orig_1, 
        dm.Phistate_sp_dest_1, 
        dm.Y_sp_dest_1,
        dm.Hstate_sp_orig_2, 
        dm.Hstate_sp_dest_2, 
        dm.Phistate_sp_orig_2, 
        dm.Phistate_sp_dest_2, 
        dm.Y_sp_dest_2,
        dm.Phi_sp_orig_1, 
        dm.Phi_sp_dest_1, 
        dm.H_sp_orig_1, 
        dm.H_sp_dest_1, 
        dm.Phi_sp_orig_2, 
        dm.Phi_sp_dest_2, 
        dm.H_sp_orig_2, 
        dm.H_sp_dest_2, 
        Insample,
        dm.ID_1, 
        dm.ID_2,
        dm.ID_sp_dest_1, 
        dm.ID_sp_dest_2,
        ID_sp_dest,
        X_orig,
        X_dest,
        MarriedHH,
        Married,
        OriginalSpouse,
        OriginalSpouse_,
        Sex_,
        s, 
        s_, 
        s_sp_orig, 
        s_sp_orig_, 
        s_sp_dest, 
        s_m_cpl, 
        s_m_cpl_, 
        s_f_cpl, 
        s_f_cpl_, 
        s_true, 
        s_sp_dest_true, 
        s_m_cpl_true, 
        s_f_cpl_true, 
        insample, 
        insample!_Bool, 
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
        h_orig_fhat,
        hstate_orig_fhat,
        x_sp_orig_fhat,
        x_sp_dest_fhat,
        xstate_sp_orig_fhat,
        xstate_sp_dest_fhat,
        y_sp_orig_fhat,
        h_sp_orig_fhat,
        hstate_sp_orig_fhat,
        x_m_cpl,
        x_f_cpl,
        xstate_m_cpl,
        xstate_m_cpl_,
        xstate_f_cpl,
        xstate_f_cpl_,
        h_m_cpl,
        h_f_cpl,
        hstate_m_cpl,
        hstate_m_cpl_,
        hstate_f_cpl,
        hstate_f_cpl_,
        x,
        xstate,
        xstate_,
        h,
        hstate,
        hstate_,
        y,
        y_m_cpl,
        y_f_cpl,
        phi,
        phistate,
        phistate_,
        phi_sp_orig,
        phistate_sp_orig,
        phistate_sp_orig_,
        h_sp_orig,
        hstate_sp_orig,
        hstate_sp_orig_,
        phi_sp_dest,
        phistate_sp_dest,
        h_sp_dest,
        hstate_sp_dest,
        id,
        firstsampled,
        firstsampled_Bool,
        I_st1,
        I_st2_single,
        I_st2_cpl,
        I_st1_,
        I_st2_single_,
        I_st2_cpl_,
        marriedhh,
        married_tminus1_,
        sex_,
        married_,
        originalspouse,
        originalspouse_reduced,
        First_t_with_sp,
        same_sp,
        NotMarried,
        NotMarried_sort,
        OriginalSpouse_sort,
        insample_sort,
        idx_NotMarried_sort,
        idx_OriginalSpouse_sort,
        idx_insample_sort,
        idx_p_st1_m_singleb4,
        idx_p_st1_f_singleb4,
        idx_p0_st1_m_singleb4,
        idx_p0_st1_f_singleb4,
        idx_p_st1_m_marriedb4,
        idx_p_st1_f_marriedb4,
        idx_p0_st1_m_marriedb4,
        idx_p0_st1_f_marriedb4,
        idx_p_st2_m_single,
        idx_p_st2_f_single,
        idx_p_st2_cpl,
        idx_p_st1_m_singleb4_reduced ,
        idx_p_st1_f_singleb4_reduced ,
        idx_p0_st1_m_singleb4_reduced ,
        idx_p0_st1_f_singleb4_reduced ,
        idx_p_st1_m_marriedb4_reduced ,
        idx_p_st1_f_marriedb4_reduced ,
        idx_p0_st1_m_marriedb4_reduced ,
        idx_p0_st1_f_marriedb4_reduced ,
        idx_p_st2_m_single_reduced ,
        idx_p_st2_f_single_reduced ,
        idx_p_st2_cpl_reduced ,
        partidx_p_st1_m_singleb4_data,
        partidx_p_st1_m_marriedb4_data,
        partidx_p_st1_f_singleb4_data,
        partidx_p_st1_f_marriedb4_data,
        partidx_p_st2_f_data,
        partidx_p_st2_m_data,
        partidx_p_st2_cpl_data,
        idx_p_st1_m_singleb4_reduced_obs,
        idx_p_st1_f_singleb4_reduced_obs,
        idx_p0_st1_m_singleb4_reduced_obs,
        idx_p0_st1_f_singleb4_reduced_obs,
        idx_p_st1_m_marriedb4_reduced_obs,
        idx_p_st1_f_marriedb4_reduced_obs,
        idx_p0_st1_m_marriedb4_reduced_obs,
        idx_p0_st1_f_marriedb4_reduced_obs,
        idx_p_st2_m_single_reduced_obs,
        idx_p_st2_f_single_reduced_obs,
        idx_p_st2_cpl_reduced_obs,
    
        idx_st1_m_singleb4_astuples,
        idx_st1_m_marriedb4_astuples,
        idx_p_st1_m_singleb4_choicesonly,
        idx_p_st1_m_marriedb4_choicesonly)

    # Also, for all these indexes, and for all arrays of integers use the smallest data type you can
    data=compress_structure(data, precision)

    return  data
end



mutable struct Data
    Adj:: Array
    Sex_1:: Array
    Phi_1:: Array
    H_1:: Array
    State_1:: Array
    Hstate_1:: Array
    Phistate_1:: Array
    Y_1:: Array
    Married_1:: Array
    Sampled_1:: Array
    Sex_2:: Array
    Phi_2:: Array 
    H_2:: Array 
    State_2:: Array 
    Hstate_2:: Array 
    Phistate_2:: Array 
    Y_2:: Array
    Married_2:: Array
    Sampled_2:: Array
    Hstate_sp_orig_1:: Array 
    Hstate_sp_dest_1:: Array 
    Phistate_sp_orig_1:: Array 
    Phistate_sp_dest_1:: Array 
    Y_sp_dest_1:: Array
    Hstate_sp_orig_2:: Array 
    Hstate_sp_dest_2:: Array 
    Phistate_sp_orig_2:: Array 
    Phistate_sp_dest_2:: Array 
    Y_sp_dest_2:: Array
    Phi_sp_orig_1:: Array 
    Phi_sp_dest_1:: Array 
    H_sp_orig_1:: Array 
    H_sp_dest_1:: Array 
    Phi_sp_orig_2:: Array 
    Phi_sp_dest_2:: Array 
    H_sp_orig_2:: Array 
    H_sp_dest_2:: Array 
    Insample:: Array
    ID_1:: Array 
    ID_2:: Array
    ID_sp_dest_1:: Array 
    ID_sp_dest_2:: Array
    ID_sp_dest:: Array
    X_orig:: Array
    X_dest:: Array
    MarriedHH:: Array
    Married:: Array
    OriginalSpouse:: Array
    OriginalSpouse_:: Array
    Sex_:: Array
    s:: Array 
    s_:: Array 
    s_sp_orig:: Array 
    s_sp_orig_:: Array 
    s_sp_dest:: Array 
    s_m_cpl:: Array 
    s_m_cpl_:: Array 
    s_f_cpl:: Array 
    s_f_cpl_:: Array 
    s_true:: Array 
    s_sp_dest_true:: Array 
    s_m_cpl_true:: Array 
    s_f_cpl_true:: Array 
    insample:: Array 
    insample!_Bool:: Array 
    td:: Array 
    t2:: Array
    t:: Array
    sex:: Array
    married:: Array
    married_tminus1:: Array
    x_orig_fhat:: Array
    x_dest_fhat:: Array
    xstate_orig_fhat:: Array
    xstate_dest_fhat:: Array
    y_orig_fhat:: Array
    h_orig_fhat:: Array
    hstate_orig_fhat:: Array
    x_sp_orig_fhat:: Array
    x_sp_dest_fhat:: Array
    xstate_sp_orig_fhat:: Array
    xstate_sp_dest_fhat:: Array
    y_sp_orig_fhat:: Array
    h_sp_orig_fhat:: Array
    hstate_sp_orig_fhat:: Array
    x_m_cpl:: Array
    x_f_cpl:: Array
    xstate_m_cpl:: Array
    xstate_m_cpl_:: Array
    xstate_f_cpl:: Array
    xstate_f_cpl_:: Array
    h_m_cpl:: Array
    h_f_cpl:: Array
    hstate_m_cpl:: Array
    hstate_m_cpl_:: Array
    hstate_f_cpl:: Array
    hstate_f_cpl_:: Array
    x:: Array
    xstate:: Array
    xstate_:: Array
    h:: Array
    hstate:: Array
    hstate_:: Array
    y:: Array
    y_m_cpl:: Array
    y_f_cpl:: Array
    phi:: Array
    phistate:: Array
    phistate_:: Array
    phi_sp_orig:: Array
    phistate_sp_orig:: Array
    phistate_sp_orig_:: Array
    h_sp_orig:: Array
    hstate_sp_orig:: Array
    hstate_sp_orig_:: Array
    phi_sp_dest:: Array
    phistate_sp_dest:: Array
    h_sp_dest:: Array
    hstate_sp_dest:: Array
    id:: Array
    firstsampled:: Array
    firstsampled_Bool:: Array
    I_st1:: Array
    I_st2_single:: Array
    I_st2_cpl:: Array
    I_st1_:: Array
    I_st2_single_:: Array
    I_st2_cpl_:: Array
    marriedhh:: Array
    married_tminus1_:: Array
    sex_:: Array
    married_:: Array
    originalspouse:: Array
    originalspouse_reduced:: Array
    First_t_with_sp:: Array
    same_sp:: Array
    NotMarried:: Array
    NotMarried_sort:: Array
    OriginalSpouse_sort:: Array
    insample_sort:: Array
    idx_NotMarried_sort:: Array
    idx_OriginalSpouse_sort:: Array
    idx_insample_sort:: Array
    idx_p_st1_m_singleb4:: Array
    idx_p_st1_f_singleb4:: Array
    idx_p0_st1_m_singleb4:: Array
    idx_p0_st1_f_singleb4:: Array
    idx_p_st1_m_marriedb4:: Array
    idx_p_st1_f_marriedb4:: Array
    idx_p0_st1_m_marriedb4:: Array
    idx_p0_st1_f_marriedb4:: Array
    idx_p_st2_m_single:: Array
    idx_p_st2_f_single:: Array
    idx_p_st2_cpl:: Array
    idx_p_st1_m_singleb4_reduced :: Array
    idx_p_st1_f_singleb4_reduced :: Array
    idx_p0_st1_m_singleb4_reduced :: Array
    idx_p0_st1_f_singleb4_reduced :: Array
    idx_p_st1_m_marriedb4_reduced :: Array
    idx_p_st1_f_marriedb4_reduced :: Array
    idx_p0_st1_m_marriedb4_reduced :: Array
    idx_p0_st1_f_marriedb4_reduced :: Array
    idx_p_st2_m_single_reduced :: Array
    idx_p_st2_f_single_reduced :: Array
    idx_p_st2_cpl_reduced :: Array
    partidx_p_st1_m_singleb4_data:: Array
    partidx_p_st1_m_marriedb4_data:: Array
    partidx_p_st1_f_singleb4_data:: Array
    partidx_p_st1_f_marriedb4_data:: Array
    partidx_p_st2_f_data:: Array
    partidx_p_st2_m_data:: Array
    partidx_p_st2_cpl_data:: Array
    idx_p_st1_m_singleb4_reduced_obs:: Array
    idx_p_st1_f_singleb4_reduced_obs:: Array
    idx_p0_st1_m_singleb4_reduced_obs:: Array
    idx_p0_st1_f_singleb4_reduced_obs:: Array
    idx_p_st1_m_marriedb4_reduced_obs:: Array
    idx_p_st1_f_marriedb4_reduced_obs:: Array
    idx_p0_st1_m_marriedb4_reduced_obs:: Array
    idx_p0_st1_f_marriedb4_reduced_obs:: Array
    idx_p_st2_m_single_reduced_obs:: Array
    idx_p_st2_f_single_reduced_obs:: Array
    idx_p_st2_cpl_reduced_obs:: Array

    idx_st1_m_singleb4_astuples:: Array
    idx_st1_m_marriedb4_astuples:: Array
    idx_p_st1_m_singleb4_choicesonly:: Array
    idx_p_st1_m_marriedb4_choicesonly:: Array
end

function comp(a,b) 
    """
    I'm just wrapping the identity function here to use "eachcol" operations... 
    """
    return a.===b 
end