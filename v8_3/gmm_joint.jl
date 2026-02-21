function gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, i, vect=0)
    # Precompute some stuff to be used in gmm_joint, for performance

    ##########################
    # CCP from logit as GMM
    ##########################
    dim_xx=size(xx)[2]
    y2_indicator=Int.(zeros((N*T2*NTypes,J)))
    for y =1:J
        y2_indicator[:,y]=(y2.==y)
    end
    # Copy xx J times horizontally
    xx_=repeat(xx',J)'
    # Copy each column of y2_indicator and U1 horizontally dim_xx times
    y2_indicator_rf=kron(y2_indicator,ones((1,dim_xx)))
    return b_all_gmm -> gmm_joint(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, xx_, y2_indicator_rf, W, vect)[i]
end

function gmm_joint(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, xx_, y2_indicator_rf, W, vect=0)
    local fvt1

    # This is a hacky way of doing this, but it's robust to passing the function to FordardDiff
    # When running a function through ForwardDiff, the type of b_all_gmm is no longer Float64 but instead it's 
    # ForwardDiff.Dual{ForwardDiff.Tag{var"#343#344", Float34}, Float34, 12}
    type=string(eltype(b_all_gmm))

    if occursin("Float64", type)
        type=Float64
        type_int = Int64
    elseif occursin("Float32", type)
        type=Float32
        type_int = Int32
    elseif occursin("Float16", type)
        type=Float16
        type_int = Int16
    end

    RX1=convert(Matrix{type}, RX1)
    Zstate=convert(Vector{type_int}, Zstate)
    Xstate=convert(Matrix{type_int}, Xstate)
    xtran=convert(Array{type}, xtran)
    xccp=convert(Array{type}, xccp)
    td=convert(Matrix{type_int}, td)
    y2=convert(Matrix{type_int}, y2)
    xx=convert(Matrix{type}, xx)
       
    # Split parameter vector. First part is for the structural equation, second part is for reduced form CCPs, third part is Pi
    b_dis=b_all_gmm[1:size(xccp)[2]+1+size(td)[2]-1]
    b1_=b_all_gmm[size(b_dis)[1]+1:size(b_dis)[1] + size(xx_)[2]]
    Π=b_all_gmm[size(b_dis)[1]+size(b1_)[1]+1:size(b_dis)[1]+size(b1_)[1]+NTypes]
    lambda=b_all_gmm[size(b_dis)[1]+size(b1_)[1]+size(Π)[1]+1]
    
    # Calculate continuation values given b1
    if vect==0
        fvt1 = fvdata(b1_, RX1, tbin, xbin, Zstate, Xstate, xtran, N, T2)
        # fvt1 = fvdata_v2(b1_, RX1, tbin, xbin, Zstate, Xstate, xtran, N, T2)
        #  fvt1, P0 = fvdata_v3(xval, zval, tbin, xbin, Zstate, Xstate, Y, PType_, xtran, N, T2)
    elseif vect==1
        fvt1 = fvdata_vect(b1_, RX1, tbin, xbin, Zstate, Xstate, xtran, N, T2)
    end

    # Calculate the q vector using the input Pi
    X_struct=hcat(xccp, reshape(fvt1,(:,1,J)), repeat(td[:, 1:T2 - 2],1,1,J))
    Like = likeCCP(b_dis, y2[index], X_struct[index,:,:])
    # Like = likeCCP_empirical(b1_, y2[index].==0, xx[index,:])
    # Like = likeCCP_empirical_nonpar(y2.==0,Xstate,Zstate, P0)
    Like2 = reshape(Like, N, T2 - 1, NTypes)
    base = reshape(prod(Like2, dims = 2),(N,NTypes))
    PType_N_2=(reshape(Π,(1,NTypes)).*base)./sum(reshape(Π,(1,NTypes)).*base,dims=2)
    # The division above can generate NaNs when working with low precision numericals
    if maximum(isnan.(PType_N_2))==true
        throw("You have NaNs in PType_N_2. Perhaps it's due to low precision")
    end
    PType_N_T_2=repeat(PType_N_2,T2)
    PType=reshape(PType_N_T_2,N*T2*NTypes)

    #######################################################################
    # Transition probabilities estimation as GMM
    #######################################################################    
    gmm_trans_prob_moment_only=gmm_trans_prob_closure(X,Z,2)
    gmm_trans_prob_momentn_only=gmm_trans_prob_closure(X,Z,3)
    moment_fhat=gmm_trans_prob_moment_only(lambda)
    res_fhat=gmm_trans_prob_momentn_only(lambda)
    

    #######################################################################
    # CCP from logit as GMM
    #######################################################################
    # Use MLE FOCs
    gmm_rf_moment_only=wlogitd_rf_gmm_closure(y2, xx, PType, W, 3)
    gmm_rf_momentn_only=wlogitd_rf_gmm_closure(y2, xx, PType, W, 4)
    moment_dis_rf=gmm_rf_moment_only(b1_)
    prod_dis_rf_weighted=gmm_rf_momentn_only(b1_)
    
    #######################################################################
    # Types probabilities Pi
    #######################################################################
    
    # To write the Pi estimation as GMM moment, I need to create (NTypes-1) vectors (one per type minus one) 
    # with one element  per observation (i.e. N*(T2-1))
    # Pretty sure this is not the correct way to do it, but it's okay for now.

    Π_=kron(reshape(Π[1:end-1],(1,NTypes-1)),ones(type_int,(N*T2,1)))
    Pi=reshape(PType,(N*T2,NTypes))
    Pi=Pi[:,1:NTypes-1]
    Pi=reshape(Pi,(N*T2,:))
    
    res_Pi=(Π_ .- Pi)
    moment_Pi=sum(res_Pi,dims=1)[:] ./(N*T2)


    #######################################################################
    # Structural equations
    #######################################################################
    if false
        # Use population type probabilities as weights
        # PType_N_T_2=reshape(PType,(N,T2,2))
        # Pi=sum(PType_N_T_2[:,1,:],dims=1)./size(PType_N_T_2,1)
        # PType_dis=repeat(Pi,NTdis)
        # PType_dis=reshape(PType_dis,NTdis*NTypes)
        # PType_dis=PType_dis[index]
        PType_dis=vcat(kron(Π_,ones(type_int,N*T2)),kron(1 .- Π_,ones(type_int,N*T2)))
    else
        # Use individual-specific type probabilities as weights
        PType_dis=PType[index]
    end

    X_dis= X_struct[index,:,:]
    Z_dis= X_dis
    Y_dis= y2[index]

    gmm_moment_only= gmm_struct_closure(Y_dis, X_dis, Z_dis, PType_dis,W,3)
    gmm_momentn_only= gmm_struct_closure(Y_dis, X_dis, Z_dis, PType_dis,W,4)
    moment_dis=gmm_moment_only(b_dis)
    prod_dis_weighted=gmm_momentn_only(b_dis)

    #######################################################################
    # Join all moments
    #######################################################################
    moment_dis_rf=moment_dis_rf[1:end-1] # Drop the moment condition for the FE on t=T2 in the reduced form
    moment=vcat(moment_dis, moment_dis_rf, moment_Pi , moment_fhat)

    # THIS IS A BIT WEIRD BECAUSE THE REDUCED FORM ESTIMATOR USED ALL DATA FROM t=1 to t=T2 WHEREAS 
    # THE STRUCTURAL ESTIMATOR ONLY USES DATA FROM t=1 TO T2-1 because you can't calculate period-ahead probabilities for T2.
    # SO the quick fix here is to just drop the data from T2 when calculating Σ
    prod_=[prod_dis_weighted prod_dis_rf_weighted[1:N*(T2-1),1:end-1] res_Pi[1:N*(T2-1),:] res_fhat[1:1:N*(T2-1)]]



    Σ=prod_' * prod_ ./ (N*(T2-1))
    # sum() only used to convert 1 by 1 matrix into scalar
    obj=sum(moment' * W * moment)

    return obj, Σ, moment, prod_
end