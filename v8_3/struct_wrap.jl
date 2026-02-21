function f_struct_wrap(bccp, dt, dt_struct_GMM, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType)

    
    bccp_st1_m_marriedb4, bccp_st1_m_singleb4, bccp_st1_f_marriedb4, bccp_st1_f_singleb4, bccp_st2_f_single, bccp_st2_m_married, bccp_st2_f_married= bccp_split(bccp, Xstruct_grid)

    bccp_dict_GMM=Dict{Int64, Union{Vector, CuArray}}()
    bccp_dict_GMM[4]=copy(bccp_st2_m_married)
    bccp_dict_GMM[3]=copy(bccp_st1_m_marriedb4)

    bccp_dict_pseudoMLE=copy(bccp_dict_GMM)
    bccp_success_dict=Dict{Int64, Optim.MultivariateOptimizationResults}()
    

    lk = ReentrantLock()
    if true
        # Estimating reduced form logit CCPs
        # This minimizes the negative pseudo-log-likelihood
        loop=[4,3]
        for i in loop
            if (i!=0)
                if bccp isa CuArray 
                    bccp_threads=CUDA.adapt(CuArray, bccp_dict_pseudoMLE[i])
                else
                    bccp_threads=bccp_dict_pseudoMLE[i]
                end
            else
                bccp_threads=bccp
            end
            # println("Start optimizing structural moment $i.")
            tic=time()
            # Using forward automatic differentiation from Optim
            # wlogit_objective_only = wlogit_rf_closure(dt, dt_rf_grid, PType, i)
            # result = optimize(wlogit_objective_only, bccp_threads, LBFGS(), o1; autodiff= :forward)

            # Using ReverseDiff/Zygote automatic differentiation
            g1 = (g,bccp_threads) -> struct_pseudoMLE_gradient_reverse!(g, bccp_threads, dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i, bccp_dict_pseudoMLE)
            wlogit_objective_only = struct_pseudoMLE_closure(dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i, bccp_dict_pseudoMLE)
            CUDA.@allowscalar result = optimize(wlogit_objective_only, g1, bccp_threads, LBFGS(), o1)
            
            
            toc=time()-tic
            # println("Done optimizing structural moment $i.")
            # println("It took $toc seconds")
            # println("")
            bccp_pseudoMLE = Optim.minimizer(result)
            # success = Optim.converged(result)
            success= result
            lock(lk) do
                bccp_dict_pseudoMLE[i]=bccp_pseudoMLE
                bccp_success_dict[i] = success
            end
        end
        # bccp_new=vcat([bccp_dict_pseudoMLE[i] for i in loop]...)
        bccp_new=vcat(bccp_dict_pseudoMLE[3] ,bccp_dict_pseudoMLE[4])
    elseif false
        for i in [4,3]
            W=I
            # Use logit approximation estimated with GMM
            if bccp isa CuArray
                bccp_threads=CUDA.adapt(CuArray, bccp_dict_GMM[i])
            else
                bccp_threads=bccp_dict_GMM[i]
            end
            println("Start optimizing reduced form moment $i.")
            tic=time()
            # wlogitd_objective_only = bccp_ -> wlogit_rf_gmm(bccp_, dt_struct_GMM, Xstruct_grid, PType, W, i, 1)[1]
            
            # Using forward differentiation from Optim (impossibly slow with the high number of params that I have)
            # println("Start optimizing moment $i with forward diff.")
            # tic=time()
            # fun=TwiceDifferentiable(wlogitd_objective_only, bccp_threads; autodiff = :forward)
            # result = optimize(fun, bccp_threads, LBFGS(), o1)
            # toc=time()-tic
            # println("Done optimizing moment $i with forward diff.")
            # println("It took $toc seconds")

            # Using Zygote or ReverseDiff
            # println("Start optimizing structural moment $i with reverse diff.")
            # tic=time()
            # g1 = (g,bccp_threads) -> gradient_struct_gmm!(g, bccp_threads ,dt_struct, Xstruct_grid, PType, W, i, 1)
            # result = optimize(wlogitd_objective_only, g1, bccp_threads, LBFGS(), o1)
            # toc=time()-tic
            # println("Done optimizing structural moment $i with reverse diff.")
            # println("It took $toc seconds")
            # println("")

            # Using ReverseDiff/Zygote automatic differentiation
            g1 = (g,bccp_threads) -> gradient_struct_gmm(g, bccp_threads, dt_struct_GMM, Xstruct_grid, fvt1_RX1, CCP_, PType, W, i)
            wlogit_objective_only = wlogit_rf_closure(dt, dt_rf_grid, dt_rf, PType, i)
            CUDA.@allowscalar result = optimize(wlogit_objective_only, g1, b1_threads, LBFGS(), o1)
            

            bccp_GMM = Optim.minimizer(result)
            success = Optim.converged(result)
            
            diff=maximum(abs.(bccp_GMM-bccp_threads))
            println(diff)
            lock(lk) do
                bccp_dict_GMM[i]=bccp_GMM
                bccp_success_dict[i] = success
            end
        end
        bccp_new=vcat([bccp_dict_GMM[i] for i in 1:5]...)
    end

    if bccp isa CuArray
        bccp_new=CUDA.adapt(CuArray, bccp_new)
    end
    
    return bccp_new, bccp_dict_pseudoMLE, bccp_success_dict
end


function gradient_struct_gmm(g, bccp_,dt_struct_GMM, Xstruct_grid, fvt1_RX1, CCP_, PType, W, i)
    wlogitd_objective_only = f_struct_gmm(bccp_, dt_struct_GMM, Xstruct_grid, fvt1_RX1, CCP_, PType, W, i)
    if bccp_ isa CuArray
        g .= Zygote.gradient(wlogitd_objective_only, bccp_)[1]
    else
        g .= ReverseDiff.gradient(wlogitd_objective_only, bccp_)  # Computes the gradient of `wlogit_objective_only` at b1
    end
end

function struct_pseudoMLE_gradient_reverse!(g, bccp_, dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i=0, bccp_dict_pseudoMLE=[])
    wlogit_objective_only = struct_pseudoMLE_closure(dt, Lambda, Xstruct_grid, fvt1_RX1, CCP_, PType, i, bccp_dict_pseudoMLE)
    if bccp_ isa CuArray
        g .= Zygote.gradient(wlogit_objective_only, bccp_)[1]
    else
        g .= ReverseDiff.gradient(wlogit_objective_only, bccp_)  # Computes the gradient of `wlogit_objective_only` at bccp
    end
end