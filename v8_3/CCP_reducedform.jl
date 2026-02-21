function CCP_reducedform(b1_, dt_rf_grid, dt_rf, i=0)
    local b1_st1_m, b1_st1_f, b1_st2_m, b1_st2_f, b1_st2_cpl, b1_split
    global dummy_size0, dummy_size10, zeros_dim_x_st2_2, zeros_dim_x_st2_cpl_2

    p0_st1_m_marriedb4=dummy_size0
    p0_st1_f_marriedb4=dummy_size0
    p0_st1_m_singleb4=dummy_size0
    p0_st1_f_singleb4=dummy_size0
    p0_st2_m=dummy_size0
    p0_st2_f=dummy_size0
    p0_st2_cpl=dummy_size0
    p_st1_m_marriedb4=dummy_size0
    p_st1_f_marriedb4=dummy_size0
    p_st1_m_singleb4=dummy_size0
    p_st1_f_singleb4=dummy_size0
    p_st2_m=dummy_size0
    p_st2_f=dummy_size0
    p_st2_cpl=dummy_size0

    # Initialize variables for storing parameter slices (as views, since iI can't allocate memory inside the optimization or Zygote will fail)
    b1_st1_m=dummy_size0 
    b1_st1_f=dummy_size0 
    b1_st2_m=dummy_size0 
    b1_st2_f=dummy_size0 
    b1_st2_cpl=dummy_size0
    
    idx1 = dt_rf.dim_x_st1_2 * dt_rf.dim_y_st1
    idx2 = dt_rf.dim_x_st2_2 * dt_rf.dim_y_st2
    idx_cpl = dt_rf.dim_x_st2_cpl_2 * dt_rf.dim_y_st2_cpl

    if i == 0
        b1_split = b1_
        b1_st1_m = get_slice(b1_split, 1, idx1)
        b1_st1_f = get_slice(b1_split, idx1+1, idx1)
        b1_st2_m = get_slice(b1_split, 2*idx1+1, idx2)
        b1_st2_f = get_slice(b1_split, 2*idx1+idx2+1, idx2)
        b1_st2_cpl = get_slice(b1_split, 2*idx1+2*idx2+1, idx_cpl)
    elseif i == 1
        b1_st1_m = b1_
    elseif i == 2
        b1_st1_f = b1_
    elseif i == 3
        b1_st2_m = b1_
    elseif i == 4
        b1_st2_f = b1_
    elseif i == 5
        b1_st2_cpl = b1_
    end

    if (i==1) | (i==0)
        # First stage parameters. 
        b1_st1_m=reshape(b1_st1_m,(dt_rf.dim_x_st1_2,phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest))

        # First stage
        Xbgrid_st1_m_marriedb4 = sum(b1_st1_m .* dt_rf_grid.RX1_st1_m_marriedb4,dims=1) ; # Dimensions 2,3,4 correspond to the chosen partner's x, h, and s. 
        Xbgrid_st1_m_marriedb4=reshape(Xbgrid_st1_m_marriedb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2)) # Now dims 1,2,3 correspond to the chosen partner, dims 4,5,6 correspond to the agent, dims 7,8,9 correspond to the previous partner

        Xbgrid_st1_m_singleb4 = sum(b1_st1_m .* dt_rf_grid.RX1_st1_m_singleb4,dims=1) # Dimensions 2,3,4 correspond to the chosen partner's x, h, and s 
        Xbgrid_st1_m_singleb4=reshape(Xbgrid_st1_m_singleb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,T2)) # Now dims 1,2,3 correspond to the chosen partner, 4,5,6 correspond to the agent

        # This is to prevent overflow
        max=maximum(Xbgrid_st1_m_marriedb4)
        expXbgrid_st1_m_marriedb4=exp.(Xbgrid_st1_m_marriedb4.-max)

        scale=1
        expXbgrid_st1_m_marriedb4_std=expXbgrid_st1_m_marriedb4 ./scale
        one_scale=exp(0 - max) / scale   # use the mean() as a scale otherwise 1/maximum() could overflow to 0
        # The baseline choice in stage 1 is always being single
        p0_st1_m_marriedb4=reshape(one_scale ./(sum(expXbgrid_st1_m_marriedb4_std, dims=(1,2,3)) .+ one_scale),(phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))
        p_st1_m_marriedb4=reshape(expXbgrid_st1_m_marriedb4_std ./(sum(expXbgrid_st1_m_marriedb4_std, dims=(1,2,3)) .+ one_scale), (phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))

        max=maximum(Xbgrid_st1_m_singleb4)
        expXbgrid_st1_m_singleb4=exp.(Xbgrid_st1_m_singleb4 .- max)
        # scale=mean(expXbgrid_st1_m_singleb4./max) .*max
        scale=1
        expXbgrid_st1_m_singleb4_std=expXbgrid_st1_m_singleb4 ./scale
        one_scale=exp(0 - max) / scale
        p0_st1_m_singleb4=reshape(one_scale ./(sum(expXbgrid_st1_m_singleb4_std, dims=(1,2,3)) .+ one_scale),(phibin_m,hbin_m,NTypes_m,T2))
        p_st1_m_singleb4=reshape(expXbgrid_st1_m_singleb4_std ./(sum(expXbgrid_st1_m_singleb4_std, dims=(1,2,3)) .+ one_scale),(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,T2))
        
    end
    if (i==2) | (i==0)
        b1_st1_f=reshape(b1_st1_f,(dt_rf.dim_x_st1_2,phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest))

        Xbgrid_st1_f_marriedb4 = sum(b1_st1_f .* dt_rf_grid.RX1_st1_f_marriedb4,dims=1) # Dimensions 2,3,4 correspond to the chosen partner's x, h, and s
        Xbgrid_st1_f_marriedb4=reshape(Xbgrid_st1_f_marriedb4,(phibin_sp_dest, hbin_sp_dest, NTypes_sp_dest, phibin_f, hbin_f, NTypes_f, phibin_sp_orig, hbin_sp_orig, NTypes_sp_orig, T2))

        Xbgrid_st1_f_singleb4 = sum(b1_st1_f .* dt_rf_grid.RX1_st1_f_singleb4,dims=1) # Dimensions 2,3,4 correspond to the chosen partner's x, h, and s
        Xbgrid_st1_f_singleb4=reshape(Xbgrid_st1_f_singleb4,(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_f,hbin_f,NTypes_f,T2))

        # This is to prevent overflow
        max=maximum(Xbgrid_st1_f_marriedb4)
        expXbgrid_st1_f_marriedb4=exp.(Xbgrid_st1_f_marriedb4 .- max)
        # scale=mean(expXbgrid_st1_f_marriedb4./max) .*max
        scale=1
        expXbgrid_st1_f_marriedb4_std=expXbgrid_st1_f_marriedb4 ./scale
        one_scale=exp(0 - max) / scale
        p0_st1_f_marriedb4=reshape(one_scale ./(sum(expXbgrid_st1_f_marriedb4_std, dims=(1,2,3)) .+ one_scale),(phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))
        p_st1_f_marriedb4=reshape(expXbgrid_st1_f_marriedb4_std ./(sum(expXbgrid_st1_f_marriedb4_std, dims=(1,2,3)) .+ one_scale),(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2))

        max=maximum(Xbgrid_st1_f_singleb4)
        expXbgrid_st1_f_singleb4=exp.(Xbgrid_st1_f_singleb4 .- max)
        # scale=mean(expXbgrid_st1_f_singleb4 ./max) .*max
        scale=1
        expXbgrid_st1_f_singleb4_std=expXbgrid_st1_f_singleb4 ./scale
        one_scale=exp(0 - max) / scale
        p0_st1_f_singleb4=reshape(one_scale ./(sum(expXbgrid_st1_f_singleb4_std, dims=(1,2,3)) .+ one_scale),(phibin_f,hbin_f,NTypes_f,T2))
        p_st1_f_singleb4=reshape(expXbgrid_st1_f_singleb4_std ./(sum(expXbgrid_st1_f_singleb4_std, dims=(1,2,3)) .+ one_scale),(phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin_m,hbin_m,NTypes_m,T2))
    end
    if (i==3) | (i==0)
        # Second stage
        # stack a vector of zeros for the baseline choice
        b1_st2_m=vcat(zeros_dim_x_st2_2,b1_st2_m)
        b1_st2_m=reshape(b1_st2_m,(dt_rf.dim_x_st2_2,J_m+1))

        # Second stage
        Xbgrid_st2_m = sum(b1_st2_m .* dt_rf_grid.RX1_st2_m,dims=1) # Dimensions 2 corresponds to the second stage choice
        Xbgrid_st2_m=reshape(Xbgrid_st2_m,(J_m+1,xbin_m,hbin_m,NTypes_m,T2)) # Now dimension 1 corresponds to the choice and dims 2,3,4 correspond to the agent's x,h, and s

        max = maximum(Xbgrid_st2_m)
        expXbgrid_st2_m=exp.(Xbgrid_st2_m .- max)
        # scale=mean(expXbgrid_st2_m ./max).*max
        scale=1
        expXbgrid_st2_m_std=expXbgrid_st2_m/scale
        one_scale=exp(0 - max) / scale

        p0_st2_m=reshape(one_scale ./(sum(expXbgrid_st2_m_std, dims=1)),(xbin_m,hbin_m,NTypes_m,T2))
        p_st2_m=reshape(expXbgrid_st2_m_std ./(sum(expXbgrid_st2_m_std, dims=1)),(J_m+1,xbin_m,hbin_m,NTypes_m,T2))
        p_st2_m=permutedims(p_st2_m,(2,3,4,1,5))
    end
    if (i==4) | (i==0)
        b1_st2_f=vcat(zeros_dim_x_st2_2,b1_st2_f)
        b1_st2_f=reshape(b1_st2_f,(dt_rf.dim_x_st2_2,J_f+1))

        Xbgrid_st2_f = sum(b1_st2_f .* dt_rf_grid.RX1_st2_f,dims=1) # Dimensions 2 corresponds to the second stage choice
        Xbgrid_st2_f=reshape(Xbgrid_st2_f,(J_f+1,xbin_f,hbin_f,NTypes_f,T2))

        max=maximum(Xbgrid_st2_f)
        expXbgrid_st2_f=exp.(Xbgrid_st2_f .- max)
        scale=1
        expXbgrid_st2_f_std=expXbgrid_st2_f/scale
        one_scale=exp(0 - max) / scale

        p0_st2_f=reshape(one_scale ./(sum(expXbgrid_st2_f_std, dims=1)),(xbin_f,hbin_f,NTypes_f,T2))
        p_st2_f=reshape(expXbgrid_st2_f_std ./(sum(expXbgrid_st2_f_std, dims=1)),(J_f+1,xbin_f,hbin_f,NTypes_f,T2))
        p_st2_f=permutedims(p_st2_f,(2,3,4,1,5))
    end
    if (i==5) | (i==0)
        b1_st2_cpl=vcat(zeros_dim_x_st2_cpl_2,b1_st2_cpl)
        b1_st2_cpl=reshape(b1_st2_cpl,(dt_rf.dim_x_st2_cpl_2,J_m+1,J_f+1))

        Xbgrid_st2_cpl = sum(b1_st2_cpl .* dt_rf_grid.RX1_st2_cpl,dims=1) # Dimensions 2 corresponds to the male's choice in the couple and dimension 3 corresponds to the female's second stage choice
        Xbgrid_st2_cpl=reshape(Xbgrid_st2_cpl,(J_m+1,J_f+1,xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))

        max=maximum(Xbgrid_st2_cpl)
        expXbgrid_st2_cpl=exp.(Xbgrid_st2_cpl .- max)
        scale=1
        expXbgrid_st2_cpl_std=expXbgrid_st2_cpl/scale
        one_scale=exp(0 - max) / scale

        p0_st2_cpl=reshape(one_scale ./(sum(expXbgrid_st2_cpl_std, dims=(1,2))),(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))
        p_st2_cpl=reshape(expXbgrid_st2_cpl_std ./(sum(expXbgrid_st2_cpl_std, dims=(1,2))),(J_m+1,J_f+1,xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))
        # Permute dimensions so that p_st2_cpl is arranged as (x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t)
        p_st2_cpl=permutedims(p_st2_cpl,(3,4,5,6,7,8,1,2,9))
    end

    #TEST: flatten all arrays
    # p0_st1_m_marriedb4=p0_st1_m_marriedb4[:]
    # p0_st1_f_marriedb4=p0_st1_f_marriedb4[:]
    # p0_st1_m_singleb4=p0_st1_m_singleb4[:]
    # p0_st1_f_singleb4=p0_st1_f_singleb4[:]
    # p0_st2_m=p0_st2_m[:]
    # p0_st2_f=p0_st2_f[:]
    # p0_st2_cpl=p0_st2_cpl[:]
    # p_st1_m_marriedb4=p_st1_m_marriedb4[:]
    # p_st1_f_marriedb4=p_st1_f_marriedb4[:]
    # p_st1_m_singleb4=p_st1_m_singleb4[:]
    # p_st1_f_singleb4=p_st1_f_singleb4[:]
    # p_st2_m=p_st2_m[:]
    # p_st2_f=p_st2_f[:]
    # p_st2_cpl=p_st2_cpl[:]

    p_RX1=P_RX1(p0_st1_m_marriedb4, p0_st1_f_marriedb4, p0_st1_m_singleb4, p0_st1_f_singleb4, p0_st2_m, p0_st2_f, p0_st2_cpl, 
    p_st1_m_marriedb4,p_st1_f_marriedb4, p_st1_m_singleb4, p_st1_f_singleb4, p_st2_m, p_st2_f, p_st2_cpl )

    return p_RX1
end

mutable struct P_RX1
    # This structure contains the CCPs for each point in the state space
    p0_st1_m_marriedb4# :: Union{Array, CuArray}
    p0_st1_f_marriedb4# :: Union{Array, CuArray}
    p0_st1_m_singleb4# :: Union{Array, CuArray}
    p0_st1_f_singleb4# :: Union{Array, CuArray}
    p0_st2_m# :: Union{Array, CuArray}
    p0_st2_f# :: Union{Array, CuArray}
    p0_st2_cpl# :: Union{Array, CuArray}
    p_st1_m_marriedb4# :: Union{Array, CuArray}
    p_st1_f_marriedb4# :: Union{Array, CuArray}
    p_st1_m_singleb4# :: Union{Array, CuArray}
    p_st1_f_singleb4# :: Union{Array, CuArray}
    p_st2_m# :: Union{Array, CuArray}
    p_st2_f# :: Union{Array, CuArray}
    p_st2_cpl# :: Union{Array, CuArray}
end

# Use views to avoid copying slices
function get_slice(data, start, len)
    return @view data[start:(start+len-1)]
end