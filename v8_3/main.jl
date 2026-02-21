using Random, LinearAlgebra,  JLD, Jacobi, MAT, SparseArrays, Distributions, CUDA, Zygote
using Optim, ForwardDiff, FiniteDiff, ReverseDiff, Optimization, OptimizationOptimJL, StatsBase
include("import.jl")

BCCP=[]
Time=[]
load_=false
MC = 1
tol = 1e-6
paramsetup()
CCP_true_sim, CCP_true=f_data_genbus4_solve()
# Generating state space for structural CCPs
Xstruct_grid=f_Xstruct_grid()
for MC in 1:50
    paramsetup()


    # Generating the data
    dt = f_data_generate_wrap(CCP_true_sim,load_)
    global N=length(dt.ID_1)
    Nobs=N*T3

    tic=time()

    # True unobserved state
    PType_true=precision.(ifelse.((dt.s_true.===dt.s).&&(dt.insample.==true),1,0).*ifelse.(((dt.s_sp_dest_true.===dt.s_sp_dest).&&(dt.married.===true)),1,ifelse.((dt.married.===false),1,0)))
    # Initialize 
    PType = rand(precision, 2*N * T2 * NTypes^2)
    PType=PType_true.+PType .* .1

    
    # True parameter vector to compare to bccp
    bccp_true=[bk... , alpha[1], alpha[2], alpha[3], alpha[4]]
    dim_bccp=Xstruct_grid.dim_Xstruct_st1_marriedb4_1+Xstruct_grid.dim_Xstruct_st2_single_1
    # Initialize parameters for structural regression
    bccp=bccp_true.+rand(precision, size(dim_bccp)) .-0.5
    
    # Initialize population type probabilities and posterior type probabilities
    Pi_cpl=reshape( precision.([0.15 0.15 ; 0.15 0.15]),(1,NTypes_m,NTypes_f))
    Pi_single_m=reshape(precision.([0.2, 0.2]),(1,NTypes_m))
    Pi_single_f=reshape(precision.([0.2, 0.2]),(1,NTypes_f))
    Pi_single=reshape(precision.([0.2, 0.2]),(1,NTypes))
    q_si_previous=zeros(precision, (N, NTypes))
    q_sij_previous=zeros(precision, (N, NTypes, NTypes))

    

    # Set to true on paramsetup if you want to run on GPU
    if RUN_ON_GPU
        # Convert arrays in structures to CuArray
        dt_CPU=dt
        dt_rf_CPU=dt_rf
        dt_rf_grid_CPU=dt_rf_grid
        Xstruct_grid_CPU=Xstruct_grid
        CCP_true_CPU=CCP_true
        dt_struct_GMM_CPU=dt_struct_GMM
        dt_CUDA, dt_rf_CUDA, dt_rf_grid_CUDA, Xstruct_grid_CUDA, CCP_true_CUDA, dt_struct_GMM_CUDA=adapt_structure_to_gpu(dt_CPU, dt_rf_CPU, dt_rf_grid_CPU, Xstruct_grid_CPU, CCP_true_CPU, dt_struct_GMM_CPU);
        dt=dt_CUDA
        dt_rf=dt_rf_CUDA
        dt_rf_grid=dt_rf_grid_CUDA
        Xstruct_grid=Xstruct_grid_CUDA
        CCP_true=CCP_true_CUDA
        dt_struct_GMM=dt_struct_GMM_CUDA

        # Convert parameters vectors to CuArray
        b1_CPU=b1
        PType_CPU=PType
        PType_true_CPU=PType_true
        xtran_CPU=xtran
        xtran_married_CPU=xtran_married
        bccp_CPU=bccp
        bccp_true_CPU=bccp_true
        Lambda_CPU=Lambda
        Pi_cpl_CPU=Pi_cpl
        Pi_single_m_CPU=Pi_single_m
        Pi_single_f_CPU=Pi_single_f

        b1=CUDA.adapt(CuArray,b1)
        PType=CUDA.adapt(CuArray,PType)
        PType_true=CUDA.adapt(CuArray,PType_true)
        xtran=CUDA.adapt(CuArray, xtran)
        xtran_married=CUDA.adapt(CuArray, xtran_married)
        bccp=CUDA.adapt(CuArray, bccp)
        bccp_true=CUDA.adapt(CuArray, bccp_true)
        Lambda=CUDA.adapt(CuArray, Lambda)
        Pi_cpl=CUDA.adapt(CuArray,Pi_cpl)
        Pi_single_m=CUDA.adapt(CuArray,Pi_single_m)
        Pi_single_f=CUDA.adapt(CuArray,Pi_single_f)
        Pi_single=CUDA.adapt(CuArray,Pi_single)
    end
    placeholders(bccp, dt)

    fvt1_RX1_true = CUDA.@allowscalar fvdata_RX1(CCP_true, xtran, xtran_married)
    fvt1_RX1=fvt1_RX1_true

    # Initialize structural CCPs
    CCP_struct=CCP_structural(bccp,Lambda, Xstruct_grid, fvt1_RX1, CCP_true)

    # Starting the EM algorithm
    # Bccp=CUDA.adapt(CuArray,zeros(precision, size(bccp')))
    Bccp=[]
    B1=[]
    condition =false
    j = 0
    W=I
    q_si_previous=CUDA.zeros(precision, (N, NTypes))
    q_sij_previous=CUDA.zeros(precision, (N, NTypes, NTypes))
    max_=1
    while !condition 

        CCP_struct=CCP_structural(bccp,Lambda, Xstruct_grid, fvt1_RX1, CCP_struct)

        Like = likeCCP_struct(dt,CCP_struct)
       
        PType_previous=PType
        Pi_cpl_previous=Pi_cpl
        Pi_single_m_previous=Pi_single_m
        Pi_single_f_previous=Pi_single_f
        Pi_single_previous=Pi_single
        
        # Update PType and Pi
        PType, q_sij, q_si=update_ptype(Like, dt, Pi_cpl, Pi_single_m, Pi_single_f)
        if false 
            # When using true PType
            PType=PType_previous
        end
        Pi_cpl, Pi_single_m, Pi_single_f=update_Pi(q_sij, q_si, dt)

        PType=precision.(PType)
        q_sij=precision.(q_sij)
        q_si=precision.(q_si)
        Pi_cpl=precision.(Pi_cpl)
        Pi_single_m=precision.(Pi_single_m)
        Pi_single_f=precision.(Pi_single_f)


        # Update the structural parameters only once in a while, to speed things up, since optimization is computationally costly
        if max_>tol*1000
            skip=10
        elseif max_>tol*100
            skip=5
        else
            skip=1
        end
        if mod(j,skip)==0
            fvt1_RX1 = CUDA.@allowscalar fvdata_RX1(CCP_struct, xtran, xtran_married)
            try
                bccp, bccp_dict_pseudoMLE, bccp_success_dict = f_struct_wrap(bccp, dt, dt_struct_GMM, Lambda, Xstruct_grid, fvt1_RX1, CCP_struct, PType)
            catch
                println("ERROR in estimation, moving on")
            end
        end

        # Stopping condition if using parametric reduced form
        if j==0
            B1=b1'
            Bccp=bccp'
        else
            B1=[B1;b1']
            maxdiff_b1=maximum(abs.(b1.-B1[j,:]))
            println("max of change in b1: $maxdiff_b1")
            Bccp=[Bccp; bccp']
            maxdiff_bccp=maximum(abs.(bccp.-Bccp[j,:]))
            println("max of change in bccp: $maxdiff_bccp")
            maxdiff_PType=maximum(abs.(PType.-PType_previous))
            println("max of change in PType: $maxdiff_PType")

            maxdiff_Pi_single_m=maximum(abs.(Pi_single_m.-Pi_single_m_previous))
            println("max of change in Pi_single_m: $maxdiff_Pi_single_m")
            maxdiff_Pi_single_f=maximum(abs.(Pi_single_f.-Pi_single_f_previous))
            println("max of change in Pi_single_f: $maxdiff_Pi_single_f")

            maxdiff_q_si=maximum(abs.(q_si.-q_si_previous))
            println("max of change in q_si: $maxdiff_q_si")

            maxdiff_q_sij=maximum(abs.(q_sij.-q_sij_previous))
            println("max of change in q_sij: $maxdiff_q_sij")

            maxdiff_Pi_cpl=maximum(abs.(Pi_cpl.-Pi_cpl_previous))
            println("max of change in Pi_cpl: $maxdiff_Pi_cpl")

            max_=max(maxdiff_bccp, maxdiff_PType, maxdiff_Pi_single_m, maxdiff_Pi_single_f, maxdiff_Pi_cpl, maxdiff_q_si, maxdiff_q_sij)
            if isnan(max_)
                println("max_ is NaN")
                break
            end
            cond1=max_<tol
            cond2=j>1000
            condition=cond1 | cond2
        end
        q_si_previous=q_si
        q_sij_previous=q_sij
        j += 1
        # println(j)
    end

    bccp, bccp_dict_pseudoMLE, bccp_success_dict = f_struct_wrap(bccp, dt, dt_struct_GMM, Lambda, Xstruct_grid, fvt1_RX1, CCP_struct, PType)
    toc=time()-tic

    println([bccp bccp_true])
    if MC==1
        BCCP=Array(bccp)'
        Time=toc
    else
        BCCP=[BCCP;Array(bccp)']
        Time=[Time,toc]
    end     
end
BCCP
    
    
