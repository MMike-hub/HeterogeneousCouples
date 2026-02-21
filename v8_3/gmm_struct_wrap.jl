function gmm_struct_wrap(binit_gmm,dt, PType,W)
    
    gmm_objective_only = gmm_struct(bccp, dt, Lambda, Xstruct_grid, fvt1_RX1, PType, W, 1)
    gmm_sigma_only= gmm_struct(bccp, dt, Lambda, Xstruct_grid, fvt1_RX1, PType, W, 2)
    gmm_moment_only= gmm_struct(bccp, dt, Lambda, Xstruct_grid, fvt1_RX1, PType, W, 3)
    gmm_momentn_only= gmm_struct(bccp, dt, Lambda, Xstruct_grid, fvt1_RX1, PType, W, 4)
    # Optimize using Optim
    fun=TwiceDifferentiable(gmm_objective_only,binit_gmm; autodiff = :forward)
    result = optimize(fun, binit_gmm, Newton(), o2)
    bccp_gmm = Optim.minimizer(result)

    # Optimize using Optimization, which allows for easier experimentation
    # fun=OptimizationFunction(gmm_objective_only, Optimization.AutoReverseDiff())
    # prob=OptimizationProblem(fun,binit_gmm)
    # result=solve(prob, BFGS())
    # bccp_gmm=result.minimizer

    D=ForwardDiff.jacobian(gmm_moment_only, bccp_gmm)
    Σ=gmm_sigma_only(bccp_gmm)
    momentn=gmm_momentn_only(bccp_gmm)
    covmat=inv(D'*W*D)*D'*W*Σ*(D'*W)'*inv(D'*W*D)'
    covmat_efficient=inv(D'*inv(Σ)*D)
    return bccp_gmm, covmat,covmat_efficient, Σ, momentn, D
end