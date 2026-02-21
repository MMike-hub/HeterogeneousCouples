function gmm_joint_wrap(bccp_gmm_1,b1, bfhat_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, Pi,W,solve=0)
    # Minimizing the vectorized version of the GMM function takes FOREVER it's unfeasible, so I do the estimation 
    # with the non-vectorized version 
    # and then only use the vectorized version to calculate the Jacobian, so that I have to run that function only once
    b_all_gmm=vcat(bccp_gmm_1,b1,Pi[:], bfhat_gmm)
    
    gmm_joint_obj_only = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 1, 0)
    gmm_joint_sigma_only = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 2, 0)
    gmm_joint_moment_only = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 3, 0)

    if solve==1
        # Solve for joint GMM estimator
        fun=TwiceDifferentiable(gmm_joint_obj_only,b_all_gmm; autodiff = :forward)
        result = optimize(fun, b_all_gmm, Newton(), o2)
        b_all_gmm = Optim.minimizer(result)
    end

    # The vectorized version if the GMM function is used to calculate the standard errors.
    # It easily goes out of memory bounds when calculating the Jacobian. Because of that, I resort to 32-bit float precision instead of 64.
    # b_all_gmm_32=convert(Vector{Float32}, b_all_gmm)

    gmm_joint_obj_only_vect = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 1, 1)
    gmm_joint_sigma_only_vect = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 2, 1)
    gmm_joint_moment_only_vect = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 3, 1)

    # Double-check that the vectorized version of the function and the loop version both have near-zero objective function 
    # value at the param estimate and near-identical Σ matrices
    obj=gmm_joint_obj_only(b_all_gmm)
    # obj_vect=gmm_joint_obj_only_vect(b_all_gmm)
    Σ=gmm_joint_sigma_only(b_all_gmm)
    # Σ_32=gmm_joint_sigma_only(b_all_gmm_32)
    # Σ_vect=gmm_joint_sigma_only_vect(b_all_gmm_32)
    if false 
        # Now calculate efficient standard errors
        W=inv(Σ) 
        gmm_joint_sigma_only = gmm_joint_closure(b_all_gmm, RX1, tbin, xbin, Zstate, Xstate, xtran, xccp, td, index, y2, xx, W, 2, 0)
        Σ=gmm_joint_sigma_only(b_all_gmm)

        D=ForwardDiff.jacobian(gmm_joint_moment_only_vect, b_all_gmm_32)
        # D=FiniteDiff.finite_difference_jacobian(gmm_joint_moment_only, b_all_gmm) 
        
        covmat=inv(D'*W*D)*D'*W*Σ*(D'*W)'*inv(D'*W*D)'
        covmat_efficient=inv(D'*inv(Matrix(Σ))*D)
    else
        covmat=[]
        covmat_efficient=[]
        D=[]
    end
    return b_all_gmm, covmat,covmat_efficient, Σ, D
end