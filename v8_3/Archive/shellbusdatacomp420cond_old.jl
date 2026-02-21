# Shell for bus program
include("genbus4.jl")
include("fvdata_RX1.jl")
include("intcond.jl")
include("intcondP.jl")
include("likeCCP.jl")
include("likeCCP_empirical.jl")
include("likeCCP_empirical_nonpar.jl")
include("plogit.jl")
include("wlogit_rf.jl")
include("wlogitd_rf.jl")
include("wlogitd_rf_gmm.jl")
include("xgrid.jl")
include("gmm_struct.jl")
include("gmm_fhat.jl")
include("reducedform_state_original.jl")
include("qupdate_logit_original.jl")
include("generate_data_wrap.jl")
include("reducedform_state_wrap.jl")
include("reducedform_choice_wrap.jl")
include("wlogit_struct.jl")
include("resCCP.jl")


using Random, LinearAlgebra,  JLD, Jacobi, MAT, SparseArrays
using Optim, ForwardDiff, FiniteDiff, ReverseDiff, Optimization, OptimizationOptimJL, StatsBase

Bccp = []
Tccp = []
Iccp = []
Binit = []
Adj = []


# Parameter values
alpha = [2, -0.15, 1, -1, 0.9, 0.4]
Beta = alpha[5]
eul = 0.577215665

# The model will be solved with a time horizon of T+10
# The first 10 periods of data will be dropped -> T2 cannot be lower than 10
# The generated dataset will go from period 11 to period T2+10, for a total panel length of T2
T = 3
T2 = T #IN THIS MODEL THERE IS NO REASON TO THROW AWAY THE FIRST 10 OBSERVATIONS. ALREADY WE HAVE RANDOM SAMPLING DOING THAT
N_ = 5
NTypes=2
NTypes_sp_orig=NTypes # ALL THESE ALTERNATIVE NAMES FOR NTypes ARE USED TO IMPROVE CODE READABILITY
NTypes_sp_dest=NTypes
NTypes_m=NTypes
NTypes_f=NTypes
NTypes_sp=NTypes
NSex=2 # Used to initiate arrays and as a memo of which dimension contains sex
NMarried = 2 # Used to initiate arrays and as a memo of which dimension contains marriage status
J=2
J_m=J # ALL THESE ALTERNATIVE NAMES FOR J ARE USED TO IMPROVE CODE READABILITY
J_f=J
J_sp_orig=J
J_sp=J


# Optimization choices
o1 = Optim.Options(store_trace = true, show_trace = true, iterations = 1000)
o2 = Optim.Options(store_trace = true, show_trace = true, iterations = 1000) #, g_tol=1e-15


# Create transition matrices
# z is a time-invariant observable characteristic
# zval = collect(0.25:0.01:1.25)
# zval = collect(0.25:0.25:1.25)
zval=[0.25; .75]
zbin = length(zval)
zbin_sp_orig=zbin # ALL THESE ALTERNATIVE NAMES FOR zbin ARE USED TO IMPROVE CODE READABILITY
zbin_sp_dest=zbin
zbin_m=zbin
zbin_f=zbin
zbin_sp=zbin
# xval = collect(0:0.25:10)
xval=Float64.(collect(0:1:5))
xbin = length(xval)
xbin_sp_orig=xbin # ALL THESE ALTERNATIVE NAMES FOR xbin ARE USED TO IMPROVE CODE READABILITY
xbin_sp_dest=xbin
xbin_m=xbin
xbin_f=xbin
xbin_sp_dest=xbin
xbin_sp_orig=xbin
xbin_orig=xbin
xbin_dest=xbin
xbin_sp=xbin
xbin_m_dest=xbin
xbin_f_dest=xbin
xbin_m_orig=xbin
xbin_f_orig=xbin
xbin_m_t1=xbin
xbin_f_t1=xbin
xbin_m_t=xbin
xbin_f_t=xbin
xbin_sp_t1=xbin
xbin_sp_t=xbin



# Size of state space
tbin = xbin * zbin

xtran = zeros(xbin, xbin, zbin, J)
xtranc = zeros(xbin, xbin, zbin, J)

# for j=1:J
#     for z in 1:zbin
#         # Notice that transitions depend on Z, observable time-invariant characteristic
#         xtran[:, :, z, j], xtranc[:, :, z, j] = xgrid(zval[z], xval , j)
#     end
# end
# xtran_1=zeros(xbin,xbin,zbin,1)
# xtranc_1=zeros(xbin,xbin,zbin,1)
# for z in 1:zbin
#     xtran_1[:,:,z,1].=repeat(xtran[1,:,z,1]',xbin)
#     xtranc_1[:,:,z,1].=repeat(xtranc[1,:,z,1]',xbin)
# end
# xtran=cat(xtran_1,xtran,dims=4)
# xtranc=cat(xtranc_1,xtranc,dims=4)
# xtran_married, xtranc_married=xgrid_married(xtran)

# The vectorized version is MUCH faster
lambda=[1.0,2.0]
xtran, xtranc, xtran_married, xtranc_married=xgrid_vect(lambda)


# z and x values for each state
zvalr = repeat(zval, inner = xbin)
xvalr = repeat(xval, outer = zbin) ./ 10


# Monte Carlos

# Starting values for FIML and CCP
alphac = [2.233, -0.1339, .4, -1, 0.9115] 
Lambda=ones(xbin_m,zbin_m,NTypes_m,xbin_f,zbin_f,NTypes_f,1,1,T2) .* 0.5 # Two singleton dimensions are for comformability to element-wise in CCP_struct 

load_=false
MC = 1
#for MC in 1:1
    # global dm, a
    # global N, T ,xtran, xtranc, xtran_married,xtranc_Married,xbin,zbin,xval,zval,T2
    # global Adj, Tccp, Bccp, Iccp, Binit, B1, b1, y2 ,xx, index, xccp, fvt1, intcondX, Pi , x , y ,z, P0, D, Σ, covmat, PType_dis
    # global  t,  x2, z2, stemp, td, PType, s2, t2, TFV, State, Zstate, Xstate , X, Z, Y, adj, base, lp, lp2, bccp, binit, PType_N_2
    # global RX1, tbin, xbin,  xtran, xtranc, Like, Like_empirical, bccp_gmm, W, X_dis, Z_dis, Y_dis, binit_gmm
    # global bccp_gmm_1, covmat_1,covmat_efficient_1, Σ_1, W_1, covmat_efficient, D_1, D, X_struct, bfhat_gmm, xtran_true, xtranc_true, bfhat_gmm

    # Generating the data
    include("genbus4.jl")
    include("generate_data_wrap.jl")
    dt = generate_data_wrap(load_)

    N=length(dt.ID_1)
    NTdis=N*T2

    tic=time()
    
    # Estimating with data CCPs
    println("Start optimization")

    
    # Transition probabilities estimation
    binit_fhat_gmm=[0.5, 0.5]
    xtran_true=copy(xtran)
    xtranc_true=copy(xtranc)
    # #Genaro's version - Returns the estimates and the transition probs.
    include("gmm_fhat.jl")
    bfhat_gmm, xtran, xtranc= gmm_fhat_wrap(binit_fhat_gmm, dt)
    #THIS ESTIMATOR SEEMS BIASED AND FOR THE LIFE OF ME I CAN'T FIGURE OUT WHY. THE FINER THE x GRID THE SMALLER THE BIAS THO

    # Setting up data for reduced form choice logit
    include("reducedform_state_wrap.jl")
    include("reducedform_state_original.jl")
    include("reducedform_state_timeinteract.jl")
    include("reducedform_state_orthocheby.jl")
    dt_rf, dt_rf_grid = reducedform_state_wrap(dt,xvalr,zvalr,zbin,xbin)
    
    # Initialize parameters for parametric reduced form estimation of choice probabilities
    b1 = 0.5 .* ones(2*dt_rf.dim_x_st1_2*dt_rf.dim_y_st1  + 2*dt_rf.dim_x_st2_2*dt_rf.dim_y_st2 +
    + dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl)
    PType = 0.5 * ones(2*N * T2 * NTypes^2)
    # PType = rand(2*N * T2 * NTypes^2)
    # PType=repeat(PType,T2)
    # PType=repeat(State,T2)
    # PType=[PType;1 .- PType]

    # To speed things up:
    if load_==true
        first_step=load("first_step.jld")
        for (key, value) in first_step
            @eval $(Symbol(key)) = $value
        end
    end

    include("wlogitd_rf_gmm.jl")
    include("reducedform_choice_wrap.jl")
    include("wlogitd_rf.jl")
    include("wlogit_rf.jl")
    include("fvdata_RX1.jl")
    include("CCP_reducedform.jl")
    include("resCCP.jl")
    b1, fvt1_RX1, CCP_rf, b1_previous = reducedform_choice_wrap(b1,dt_rf,PType)
  
    # Starting the EM algorithm
    j = 0
    

    include("Xstruct.jl") 
    Xstruct_grid=Xstruct()
    
    bccp = vcat(alphac[:], zeros(T2 - 2))
    dim_bccp=2*Xstruct_grid.dim_Xstruct_st1_marriedb4_1+2*Xstruct_grid.dim_Xstruct_st1_singleb4_1+2*Xstruct_grid.dim_Xstruct_st2_single_1+2*Xstruct_grid.dim_Xstruct_st2_married_1
    bccp=ones(dim_bccp)
    
    binits = zeros(3)
    binit = binits
    index = dt.t2 .< 1
    lp = []
    lp2= []
    if load_==true
        cond = 1 
    else
        cond = 0
    end
    tol = 1e-7
    W=I


    B1=[]
    
    while cond == 0
        # Updating PType
        include("resCCP.jl")
        include("CCP_struct.jl")
        CCP_struct=CCP_structural(bccp,Lambda, Xstruct_grid, fvt1_RX1, CCP_rf)
        include("likeCCP.jl")
        Like = likeCCP(dt,CCP_struct)
        # Like = likeCCP_empirical(b1, y2[index].==0, xx[index,:])
        # Like = likeCCP_empirical_nonpar(y2.==0,Xstate,Zstate, P0)
        if false #& (j % 10 == 0)
            include("likeCCP_empirical.jl")
            Like_empirical = likeCCP_empirical(b1, y2[index], xx[index,:])
            include("likeCCP_empirical_nonpar.jl")
            Like_empirical = likeCCP_empirical_nonpar(y2[index],Xstate[:,1:T2-1],Zstate, P0)
            # println("s=1")
            # println((Like_empirical-Like)[1:5])
            # println("s=2")
            # println((Like_empirical-Like)[NT2+1:NT2+5])
            println(maximum(abs.(Like_empirical-Like)))
        end
        if false 
            Like2 = reshape(Like, N, T2 - 1, NTypes)
            base = reshape(prod(Like2, dims = 2),(N,NTypes))

            # This is the fully nonparametric Q update
            PType_previous=PType
            PType_N_T_2=reshape(PType,(N,T2,2))
            Pi=sum(PType_N_T_2[:,1,:],dims=1)./size(PType_N_T_2,1)
            PType_N_2=(Pi.*base)./sum(Pi.*base,dims=2)
            PType_N_T_2=repeat(PType_N_2,T2)
            PType=reshape(PType_N_T_2,N*T2*2)

            lp=vcat(lp,sum((PType.-PType_previous).^2))
        else
            include("update_ptype.jl")
            ID_sp = vcat(dt.ID_sp_dest_1,dt.ID_sp_dest_2)
            ID_sp = get.(coalesce.(ID_sp,Ref((0,0))), 1, missing)
            PType = update_ptype(Like, ID_sp, dt.Insample) 
        end

        b1, fvt1_RX1, CCP_rf, b1_previous = reducedform_choice_wrap(b1,dt_rf,PType)

        # NEED TO USE THIS AS LONG AS THE Like ARRAY IS CALCULATED USING likeCCP
        # Original MLE estimator        
        include("wlogit_struct.jl")
        wlogit_objective = wlogit_struct_closure(dt, Lambda, Xstruct_grid, fvt1_RX1, PType)
        fun=TwiceDifferentiable(wlogit_objective,bccp; autodiff = :forward)
        result = optimize(fun, bccp, BFGS(), o1)
        bccp = Optim.minimizer(result)

        if true
            # Stopping condition if using parametric reduced form
            if j==0
                B1=b1[1:4]'
            else
                B1=[B1;b1[1:4]']
                println(b1[1:4].-B1[j,:])
                lp2=sum((b1[1:4].-B1[j,:]).^2)
                junk=lp[j+1]<tol
                junk2=lp2<tol
                cond=junk*junk2+(j>150)
            end
        elseif false
            # Stopping condition if using nonparametric reduced form
            if j==0
                lp=maximum(P0.-P0_previous)
            else
                lp=vcat(lp, maximum(P0.-P0_previous))
            end
            cond=lp[end]<tol
        end
            
      

        j += 1
        println(j)
    end
    
    save("first_step.jld", "b1", b1, "PType", PType, "Pi", Pi)
    save("busdata0210_dataset.jld", "Zstate", Zstate, "Xstate", Xstate, "Y",  Y , "X", X, "Z", Z, "y2", y2, "x", x, "x2",  x2, "z", z , "z2", z2 , "s2",  s2, "td",  td, "t2",  t2) 

    if true
        # Estimate structural parameters via MLE
        wlogit_objective = wlogit_struct_closure(dt, Lambda, Xstruct_grid, fvt1_RX1, PType)
        fun=TwiceDifferentiable(wlogit_objective,bccp; autodiff = :forward)
        result = optimize(fun, bccp, BFGS(), o1)
        bccp = Optim.minimizer(result)

    end

    if true
        # Estimate structural parameters via GMM
        CCP_struct=CCP_structural(bccp,Lambda, Xstruct_grid, fvt1_RX1)
        # fvt1 = fvdata_v2(b1, RX1, tbin, xbin, Zstate, Xstate, xtran, N, T2)
        #  fvt1, P0 = fvdata_v3(xval, zval, tbin, xbin, Zstate, Xstate, Y, PType, xtran, N, T2)
        
        # include("gmm_struct.jl")
        include("gmm_struct_v2.jl")
        include("gmm_struct_wrap.jl")
        include("wlogit_struct.jl")
        include("wlogitd_rf_gmm.jl")
        include("wlogitd_rf.jl")
        W_1=I
        bccp_gmm_1, covmat_1,covmat_efficient_1, Σ_1, momentn_1, D_1 =gmm_struct_wrap(binit_gmm,dt, PType,W)
        W_1=inv(Σ_1)
        bccp_gmm_1, covmat_1,covmat_efficient_1, Σ_1, momentn_1, D_1=gmm_struct_wrap(binit_gmm,dt, PType,W)

        # This adjusts standard error of structural parameters for firs-step incidental parameters without running joint GMM estimator
        # and just taking the parameters as given from first and second step
        include("gmm_joint_wrap.jl")
        include("wlogitd_rf_gmm.jl")
        include("gmm_joint.jl")
        include("gmm_struct.jl")
        include("fvdata_RX1.jl")
        include("fvdata_dt.jl")
        W=I
        bccp_gmm, covmat,covmat_efficient, Σ, D=gmm_joint_wrap(bccp_gmm_1,b1, bfhat_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, Pi,W,0)
        # first_step=load("x1_idx.jld")
        # for (key, value) in first_step
        #     @eval $(Symbol(key)) = $value
        # end
        
    end


    tccp = time() - tic

    if MC==1
        Bccp=bccp'
        Binit= binit'
        # For some reason, stacking Adj doesn't work. It's identical to all other vectors I'm stacking, but it doesn't work
        # Adj = adj'
        Tccp=tccp
        Iccp = j
    else
        Bccp = vcat(Bccp, bccp')
        Binit = vcat(Binit, binit')
        # Adj = vcat(Adj, adj')
        Tccp = vcat(Tccp, tccp)
        Iccp = vcat(Iccp, j)
    end

    b_true=vcat(alpha,Adj[2:T2-1]) 

    println("MC iteration $MC complete")
    println(bccp[1:4])
    println(bccp_gmm[1:4])
    # save("busdata0210.jld", "Bccp", Bccp, "Tccp", Tccp, "Iccp", Iccp, "Binit", Binit, "Adj", Adj)
#end
