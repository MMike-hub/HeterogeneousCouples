# THE REDUCED FORM LOGIT IS CHARACTERIZED BY CHOICE-SPECIFIC PARAMETER VALUES
# unlike the structural logit, where the parameter values do not necessarily change across choices

# Remember that this MLE combines many distinct types of CCPs, 
# 1) First-stage CCPs, one for males and one for females
# 2) Second-stage CCPs, one for single males, one for single females, and one for couples.

struct Precompute_reducedform
    idx_st1_m:: Vector{CartesianIndex{4}}
    idx_st1_f:: Vector{CartesianIndex{4}}
    idx_st2_single_m:: Vector{CartesianIndex{2}}
    idx_st2_single_f:: Vector{CartesianIndex{2}}
    idx_st2_married:: Vector{CartesianIndex{3}}
end

function wlogit_rf_closure(dt_rf, PType)
    # For performance, I took the definition of idx out of the wlogit function

    xstate_sp_dest_=coalesce.(dt.xstate_sp_dest[:],1)
    zstate_sp_dest_=coalesce.(dt.zstate_sp_dest[:],1)
    y_=coalesce.(dt_rf.y[:],1)
    y_m_cpl_=coalesce.(dt_rf.y_m_cpl[:],1)
    y_f_cpl_=coalesce.(dt_rf.y_f_cpl[:],1)



    #First stage
    idx_st1_m = CartesianIndex.(collect(1:dt_rf.dim_x_st1_1), xstate_sp_dest_, zstate_sp_dest_, dt.s_sp_dest[:])
    idx_st1_f = CartesianIndex.(collect(1:dt_rf.dim_x_st1_1), xstate_sp_dest_, zstate_sp_dest_, dt.s_sp_dest[:])


    #Second stage
    idx_st2_single_m = CartesianIndex.(collect(1:dt_rf.dim_x_st2_1), y_)
    idx_st2_single_f = CartesianIndex.(collect(1:dt_rf.dim_x_st2_1), y_)
    idx_st2_married = CartesianIndex.(collect(1:dt_rf.dim_x_st2_cpl_1), y_m_cpl_,  y_f_cpl_)

    prec_rf=Precompute_reducedform(idx_st1_m,idx_st1_f,idx_st2_single_m,idx_st2_single_f,idx_st2_married)

    println("IDX evaluated")
    return b1 -> wlogit_rf(b1,  dt_rf, prec_rf, PType, J)
end

function wlogit_rf(b1,  dt_rf, prec_rf, PType, J)


    #For testing
    # b1=ones(dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl)
    
    b1_=b1

    b1_st1_m=b1_[1:dt_rf.dim_x_st1_2*dt_rf.dim_y_st1]
    b1_=b1_[dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+1:end]
    b1_st1_f=b1_[1:dt_rf.dim_x_st1_2*dt_rf.dim_y_st1]
    b1_=b1_[dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+1:end]
    b1_st2_m=b1_[1:dt_rf.dim_x_st2_2*dt_rf.dim_y_st2]
    b1_=b1_[dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+1:end]
    b1_st2_f=b1_[1:dt_rf.dim_x_st2_2*dt_rf.dim_y_st2]
    b1_=b1_[dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+1:end]
    b1_st2_cpl=b1_[1:dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl]
    b1_=b1_[dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl+1:end] #if b1_ should be empty at this point

    xx_st1_=coalesce.(dt_rf.xx_st1,1)
    xx_st2_=coalesce.(dt_rf.xx_st2,1)
    xx_st2_cpl_=coalesce.(dt_rf.xx_st2_cpl,1)
    married_=coalesce.(dt_rf.married,1)
    sex_=coalesce.(dt_rf.sex,1)


    # First stage
    b1_st1_m=reshape(b1_st1_m,(dt_rf.dim_x_st1_2,1,xbin,zbin,NTypes))
    b1_st1_f=reshape(b1_st1_f,(dt_rf.dim_x_st1_2,1,xbin,zbin,NTypes))
    # these Xb_st1 DO NOT contain the baseline choice (becoming/remaining single)
    Xb_st1_m = sum(xx_st1_' .* b1_st1_m,dims=1)
    Xb_st1_m=reshape(Xb_st1_m,(:,xbin,zbin,NTypes))

    Xb_st1_f = sum(xx_st1_' .* b1_st1_f,dims=1)
    Xb_st1_f=reshape(Xb_st1_f,(:,xbin,zbin,NTypes))

    # Second stage
    # stack a vector of zeros for the baseline choice
    b1_st2_m=vcat(zeros(dt_rf.dim_x_st2_2),b1_st2_m)
    b1_st2_m=reshape(b1_st2_m,(dt_rf.dim_x_st2_2,1,J+1))

    b1_st2_f=vcat(zeros(dt_rf.dim_x_st2_2),b1_st2_f)
    b1_st2_f=reshape(b1_st2_f,(dt_rf.dim_x_st2_2,1,J+1))

    b1_st2_cpl=vcat(zeros(dt_rf.dim_x_st2_cpl_2),b1_st2_cpl)
    b1_st2_cpl=reshape(b1_st2_cpl,(dt_rf.dim_x_st2_cpl_2,1,J+1,J+1))

    # All tese Xb_st2 arrays contain the the baseline choice too along dimension 3
    Xb_st2_m = sum(xx_st2_' .* b1_st2_m,dims=1)
    Xb_st2_m = reshape(Xb_st2_m,(:,J+1))
    Xb_st2_f = sum(xx_st2_' .* b1_st2_f,dims=1)
    Xb_st2_f = reshape(Xb_st2_f,(:,J+1))
    Xb_st2_cpl = sum(xx_st2_cpl_' .* b1_st2_cpl,dims=1)
    Xb_st2_cpl = reshape(Xb_st2_cpl,(:,J+1,J+1))



    
    # The likelihood is the product of first-stage and second-stage for each individual at each period
    Like = PType' * ((dt_rf.insample).*(
    (log.(sum(exp.(Xb_st1_m),dims=(2,3,4)).+1) .- Xb_st1_m[prec_rf.idx_st1_m].*married_ .-(married_.==0)).*(sex_.==1) .+
    (log.(sum(exp.(Xb_st1_f),dims=(2,3,4)).+1) .- Xb_st1_f[prec_rf.idx_st1_f].*married_ .-(married_.==0)).*(sex_.==2) .+
    (log.(sum(exp.(Xb_st2_m),dims=2)) .- Xb_st2_m[prec_rf.idx_st2_single_m]).*(sex_.==1).*(married_.==0) .+
    (log.(sum(exp.(Xb_st2_f),dims=2)) .- Xb_st2_f[prec_rf.idx_st2_single_f]).*(sex_.==2).*(married_.==0) .+
    (log.(sum(exp.(Xb_st2_cpl),dims=(2,3))) .- Xb_st2_cpl[prec_rf.idx_st2_married]).*(married_.==1))
    .*(dt_rf.insample.==1))[:]

    # This is just to ensure that Like is a scalar
    Like=sum(Like)
    return Like
end