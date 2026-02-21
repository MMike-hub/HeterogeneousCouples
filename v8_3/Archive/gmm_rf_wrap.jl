function gmm_rf_wrap(b1_previous, y2, xx, PType)
    wlogitd_objective = b -> wlogitd_gmm(b,y2 .== 0, xx, PType,W)
    wlogitd_objective_only = b -> wlogitd_objective(b)[1];
    wlogitd_sigma_only = b -> wlogitd_objective(b)[2]; 
    wlogitd_moment_only = b -> wlogitd_objective(b)[3]; 
    wlogitd_momentn_only = b -> wlogitd_objective(b)[4]; 
    fun=TwiceDifferentiable(wlogitd_objective_only,b1_previous; autodiff = :forward)
    result = optimize(fun, b1_previous, BFGS(), o2)
    b1_GMM = Optim.minimizer(result)
    Σ_rf=wlogitd_sigma_only(b1_GMM)
    momentn_rf=wlogitd_momentn_only(b1_GMM)
    D_rf=ForwardDiff.jacobian(gmm_moment_only, b1_GMM)
    return b1_GMM, Σ_rf, momentn_rf, D_rf
end