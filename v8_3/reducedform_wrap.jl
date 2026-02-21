function f_reducedform_wrap(b1,dt_rf_grid, dt_rf, PType)
    local b1_st1_m, b1_st1_f, b1_st2_m, b1_st2_f, b1_st2_cpl
    global b1_success_dict, b1_previous_inner

    b1_previous=b1

    b1_split=b1

    b1_st1_m=b1_split[1:dt_rf.dim_x_st1_2*dt_rf.dim_y_st1] # Parameters concerning 1st-stage choice of male
    b1_split=b1_split[dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+1:end] 
    b1_st1_f=b1_split[1:dt_rf.dim_x_st1_2*dt_rf.dim_y_st1] # Parameters concerning 1st-stage choice of female
    b1_split=b1_split[dt_rf.dim_x_st1_2*dt_rf.dim_y_st1+1:end]
    b1_st2_m=b1_split[1:dt_rf.dim_x_st2_2*dt_rf.dim_y_st2] # Parameters concerning 2nd-stage choice of male
    b1_split=b1_split[dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+1:end]
    b1_st2_f=b1_split[1:dt_rf.dim_x_st2_2*dt_rf.dim_y_st2] # Parameters concerning 2nd-stage choice of female
    b1_split=b1_split[dt_rf.dim_x_st2_2*dt_rf.dim_y_st2+1:end]
    b1_st2_cpl=b1_split[1:dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl] # Parameters concerning 2nd-stage choice of couple
    b1_split=b1_split[dt_rf.dim_x_st2_cpl_2*dt_rf.dim_y_st2_cpl+1:end] # b1_split should be empty at this point

    b1_dict_GMM=Dict{Int64, Union{Vector, CuArray}}()
    b1_dict_GMM[1]=b1_st1_m
    b1_dict_GMM[2]=b1_st1_f
    b1_dict_GMM[3]=b1_st2_m
    b1_dict_GMM[4]=b1_st2_f
    b1_dict_GMM[5]=b1_st2_cpl

    b1_dict_pseudoMLE=copy(b1_dict_GMM)
    b1_success_dict=Dict{Int64, Optim.MultivariateOptimizationResults}()
    
    

    lk = ReentrantLock()
    if true
        # Estimating reduced form logit CCPs
        # This minimizes the negative presudo-log likelihood
        for i in 1:5
            if b1 isa CuArray
                b1_threads=CUDA.adapt(CuArray, b1_dict_pseudoMLE[i])
            else
                b1_threads=b1_dict_pseudoMLE[i]
            end
            println("Start optimizing reduced form moment $i with reverse diff.")
            tic=time()
            # Using forward automatic differentiation from Optim
            # wlogit_objective_only = rf_pseudoMLE_closure(dt, dt_rf_grid, PType, i)
            # result = optimize(wlogit_objective_only, b1_threads, LBFGS(), o1; autodiff= :forward)

            # Using ReverseDiff/Zygote automatic differentiation
            g1 = (g,b1_threads) -> rf_pseudoMLE_gradient_reverse!(g, b1_threads, dt, dt_rf_grid, dt_rf, PType, i)
            wlogit_objective_only = rf_pseudoMLE_closure(dt, dt_rf_grid, dt_rf, PType, i)
            CUDA.@allowscalar result = optimize(wlogit_objective_only, g1, b1_threads, GradientDescent(), o1)
            
            # Using hard-coded gradient. THIS IS BUGGY. I PROBABLY MADE SOME MISTAKES IN CODING THE GRADIENT
            # fg! = rf_pseudoMLE_closure(dt, dt_rf_grid, dt_rf, PType, i);
            # b1_previous_inner=b1_threads
            # # CUDA.@allowscalar result = Optim.optimize(Optim.only_fg!(fg!), b1_threads, Optim.LBFGS(linesearch = Optim.BackTracking()), o1)
            # CUDA.@allowscalar result = Optim.optimize(Optim.only_fg!(fg!), b1_threads, Optim.GradientDescent(), o1)
            toc=time()-tic
            println("Done optimizing reduced form moment $i with reverse diff.")
            println("It took $toc seconds")
            println("")
            b1_pseudoMLE = Optim.minimizer(result)
            # success = Optim.converged(result)
            success= result
            lock(lk) do
                b1_dict_pseudoMLE[i]=b1_pseudoMLE
                b1_success_dict[i] = success
            end
        end
        b1_new=vcat([b1_dict_pseudoMLE[i] for i in 1:5]...)
    elseif false
        for i in 1:5
            W=I
            # Use logit approximation estimated with GMM
            if b1 isa CuArray
                b1_threads=CUDA.adapt(CuArray, b1_dict_GMM[i])
            else
                b1_threads=b1_dict_GMM[i]
            end
            println("Start optimizing reduced form moment $i.")
            tic=time()
            # wlogitd_objective_only = b1_threads -> wlogit_rf_gmm(b1_threads, dt_rf, dt_rf_grid, PType, W, i, 1)[1]
            
            # Using forward differentiation from Optim (impossibly slow with the high number of params that I have)
            # println("Start optimizing moment $i with forward diff.")
            # tic=time()
            # fun=TwiceDifferentiable(wlogitd_objective_only, b1_threads; autodiff = :forward)
            # result = optimize(fun, b1_threads, LBFGS(), o1)
            # toc=time()-tic
            # println("Done optimizing moment $i with forward diff.")
            # println("It took $toc seconds")

            # Using explicit gradient function
            # println("Start optimizing reduced form moment $i with reverse diff.")
            # tic=time()
            # g1 = (g,b1_threads) -> rf_gmm_gradient_reverse!(g, b1_threads ,dt_rf, dt_rf_grid, PType, W, i, 1)
            # result = optimize(wlogitd_objective_only, g1, b1_threads, LBFGS(), o1)
            # toc=time()-tic
            # println("Done optimizing reduced form moment $i with reverse diff.")
            # println("It took $toc seconds")
            # println("")

            # Using ReverseDiff/Zygote automatic differentiation
            g1 = (g,b1_threads) -> rf_gmm_gradient_reverse!(g, b1_threads ,dt_rf, dt_rf_grid, PType, W, i)
            wlogit_objective_only = b1_ -> wlogit_rf_gmm(b1_, dt_rf, dt_rf_grid, PType, W, i, 1)
            CUDA.@allowscalar result = optimize(wlogit_objective_only, g1, b1_threads, LBFGS(), o1)
            

            toc=time()-tic
            println("Done optimizing reduced form moment $i.")
            println("It took $toc seconds")
            println("")
            b1_GMM = Optim.minimizer(result)
            # success = Optim.converged(result)
            success=result
            
            diff=maximum(abs.(b1_GMM-b1_threads))
            println(diff)
            lock(lk) do
                b1_dict_GMM[i]=b1_GMM
                b1_success_dict[i] = success
            end
        end
        b1_new=vcat([b1_dict_GMM[i] for i in 1:5]...)
    end

    if b1 isa CuArray
        b1_new=CUDA.adapt(CuArray, b1_new)
    end

    # Calculating FV terms
    CCP_rf=CCP_reducedform(b1, dt_rf_grid, dt_rf)
    fvt1_RX1 = CUDA.@allowscalar fvdata_RX1(CCP_rf, xtran, xtran_married)

    return b1_new, fvt1_RX1, CCP_rf, b1_previous, b1_success_dict
end

function rf_gmm_gradient_reverse!(g, b1_,dt_rf, dt_rf_grid, PType, W, i=0)
    wlogitd_objective_only = b1_ -> reducedform_gmm(b1_, dt_rf, dt_rf_grid, PType, W, i, 1)
    if b1_ isa CuArray
        g .= Zygote.gradient(wlogitd_objective_only, b1_)[1]
    else
        g .= ReverseDiff.gradient(wlogitd_objective_only, b1_)  # Computes the gradient of `wlogit_objective_only` at b1
    end
end

function rf_pseudoMLE_gradient_reverse!(g,b1_, dt, dt_rf_grid, dt_rf, PType, i=0)
    wlogit_objective_only = rf_pseudoMLE_closure(dt, dt_rf_grid, dt_rf, PType, i)
    if b1_ isa CuArray
        g .= Zygote.gradient(wlogit_objective_only, b1_)[1]
    else
        g .= ReverseDiff.gradient(wlogit_objective_only, b1_)  # Computes the gradient of `wlogit_objective_only` at b1
    end
end