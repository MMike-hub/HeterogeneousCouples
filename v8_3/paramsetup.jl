# Parameter values

function paramsetup()
    global alpha, Beta, eul, T, T2, N_, NTypes, NTypes_sp_orig, NTypes_sp_dest, NTypes_m, NTypes_f, NTypes_sp, 
    NSex, NMarried, J, J_m, J_f, J_sp_orig, J_sp, o1, o2, hval, hbin, hbin_sp_orig, hbin_sp_dest, hbin_m, hbin_f, hbin_m_orig, hbin_f_orig,
    hbin_sp, xval, zval, phival, xbin, zbin, phibin, xstate, zstate, phistate_list, phibin_sp_orig, phibin_sp_dest, phibin_m, phibin_f, 
    phibin_orig, phibin_dest, phibin_sp, phibin_m_dest, phibin_f_dest, phibin_m_orig, phibin_f_orig, 
    phibin_m_t1, phibin_f_t1, phibin_m_t, phibin_f_t, phibin_sp_t1, phibin_sp_t, xtran, xtranc, xtran_married, 
    xtranc_married, alphac, Lambda, Lambda_T, xbin_sp_orig,xbin_sp_dest, xbin_m, xbin_f, xbin_sp_dest, xbin_sp_orig, 
    xbin_orig, xbin_dest, xbin_sp, xbin_m_dest, xbin_f_dest, xbin_m_orig, xbin_f_orig, xbin_m_t1, xbin_f_t1, xbin_m_t, xbin_f_t, 
    xbin_sp_t1, xbin_sp_t, bk, phival_grid, phistate_grid, phidim, xdim, zdim, bpkshift,
    T3, P_drop, P_sample, COND_ORIG_SPOUSE, precision, RUN_ON_GPU, FLEX_EST, P_male, Adj, Pi_m, Pi_f, eul, Beta

    precision=Float64

    T = 20
    T2 = 10 #IN THIS MODEL THERE IS NO REASON TO THROW AWAY THE FIRST 10 OBSERVATIONS. ALREADY WE HAVE RANDOM SAMPLING DOING THAT
    T3=10
    N_ = 5000
    NTypes=2

    RUN_ON_GPU= true
    P_drop=0.0
    P_sample=1.0
    P_male=0.5
    eul = precision.(0.577215665)
    Beta = precision.(0.6)
    Pi_m=0.4
    Pi_f=0.4

    alpha = precision.([2, -0.15, 1, -1, Beta, Pi_m, Pi_f]) #[2, -0.15, 1, -1, 0.0, 0.4] #
    bk = precision.([ -4] )
    
    
    # Adj is a set of year fixed effects
    Adj = zeros(precision,T)
    # t = 2
    # while t <= T 
    #     Adj[t] = 0.7 * Adj[t-1] + 0.5 * randn()
    #     t += 1
    # end

    COND_ORIG_SPOUSE=false #condition on origin spouse
    # If FLEX_EST==true then the structural parameters for 1st-stage men, 1st-stage women, 2nd-stage men, 2nd-stage women are allowed to be all different
    # If FLEX_EST==false then the structural parameter are restricted to match the DGP (i.e. same utility parameters for men and women)
    FLEX_EST=false


     
    NTypes_sp_dest=NTypes # ALL THESE ALTERNATIVE NAMES FOR NTypes ARE USED TO IMPROVE CODE READABILITY
    
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
    o1 = Optim.Options(store_trace = true, show_trace = false, iterations = 1000, g_tol=1e-10)
    o2 = Optim.Options(store_trace = true, show_trace = true, iterations = 1000) #, g_tol=1e-15


    # Create transition matrices
    # h is a time-invariant observable characteristic
    # hval = precision.(collect(0.25:0.01:1.25))
    # hval = precision.(collect(0.25:0.25:1.25))
    # hval=precision.([0.25; .75])
    hval=precision.([1.0])
    hbin = length(hval)
    # ALL THESE ALTERNATIVE NAMES FOR hbin ARE USED TO IMPROVE CODE READABILITY
    hbin_sp_dest=hbin
    hbin_m=hbin
    hbin_f=hbin
    hbin_sp=hbin
    hbin_f_orig=hbin
    hbin_m_orig=hbin

    # x is a set of time-variant observable characteristics
    # xval=precision.(collect(0:.125:25))
    xval=precision.(collect(0:1:25))
    # xval=precision.(collect(0:1:2))
    
    # zval=precision.(collect(0:.5:2.5))
    zval=precision.(collect(0:1))
    # zval=precision.([0.0])
    xbin=length(xval)
    zbin=length(zval)
    xstate=collect(1:xbin)
    zstate=collect(1:zbin)
    phival_grid= Array{Matrix{precision}}(undef, xbin, zbin)
    phistate_grid=Array{Matrix{Int64}}(undef, xbin, zbin)

    for i in 1:xbin
        for j in 1:zbin
            phival_grid[i, j] = [xval[i] zval[j]]
            phistate_grid[i, j] = [i j]
        end
    end
    phidim=length(phival_grid[1,1])
    xdim=1
    zdim=1
    phival=vcat(phival_grid[:]...)
    phistate_list=vcat(phistate_grid[:]...)
    phibin = size(phival,1)
    # ALL THESE ALTERNATIVE NAMES FOR phibin ARE USED TO IMPROVE CODE READABILITY
    phibin_sp_dest=phibin
    phibin_m=phibin
    phibin_f=phibin
    phibin_sp_dest=phibin
    phibin_orig=phibin
    phibin_dest=phibin
    phibin_sp=phibin
    phibin_m_dest=phibin
    phibin_f_dest=phibin
    phibin_m_orig=phibin
    phibin_f_orig=phibin
    phibin_m_t1=phibin
    phibin_f_t1=phibin
    phibin_m_t=phibin
    phibin_f_t=phibin
    phibin_sp_t1=phibin
    phibin_sp_t=phibin

     # ALL THESE ALTERNATIVE NAMES FOR phibin ARE USED TO IMPROVE CODE READABILITY
    xbin_sp_dest=xbin
    xbin_m=xbin
    xbin_f=xbin
    xbin_sp_dest=xbin
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

    # If we use a model where the previous spouse's characteristics don't matter in first-stage choices, that's equivalent to 
    # pretending that in the first stage, previous spouses are always in any single state. So I set all state variables to 1 
    if COND_ORIG_SPOUSE==0
        NTypes_sp_orig=1
        xbin_sp_orig=1
        phibin_sp_orig=1
        hbin_sp_orig=1 
    else
        NTypes_sp_orig=NTypes
        xbin_sp_orig=xbin
        phibin_sp_orig=phibin
        hbin_sp_orig=hbin
    end
    
    b_fhat=precision.([1.0,2.0])
    bpkshift=precision.([1.0 2.0])
    xtran, xtranc, xtran_married, xtranc_married=xgrid_vect(b_fhat, bpkshift)


    # Starting values for FIML and CCP
    alphac = [2.233, -0.1339, .4, -1, 0.9115] 
    Lambda=precision.(ones(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,1,1,T2) .* 0.5) # Two singleton dimensions are for comformability to element-wise in genbus4 
    Lambda_T=precision.(ones( xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,1,1,T) .* 0.5)
end