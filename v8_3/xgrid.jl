function xgrid(b_fhat, bpkshift,  phival, j)
    xub = vcat(xval[2:end], Inf)
    xtran = zeros(precision, xbin_orig, phibin_dest, hbin, J)
    xtranc = zeros(precision, xbin_orig, phibin_dest, hbin, J)
    
    
    for j in 1:J
        for h in 1:hbin
            lcdf2= 0.0
            for z_dest in 1:zbin
                lcdf1 = 0.0
                for x_dest in 1:xbin
                    # Consider the different discrete choices as partial engine replacements / maintenance.
                    # The larger the value of j, the lower the probability that the next period mileage is high
                    # i.e. the larger the j, the more thorough the maintenance / engine parts replacement
                    p_mileage_transition=(xub[x_dest] .>= xval) .* (1 .- exp.(-hval[h] .* b_fhat[j] .* (xub[x_dest] .- xval)) .- lcdf1)
                    # Probability of instrument conditional on mileage. Assume it's binomial with zbin-1 trials (because in the binomial 0 is a possible result) and p=1/(exp(bpkshift[1]*x+bpkshift[2]*hval[h])+1)
                    p_= 1 ./(exp.(xval[x_dest] .* bpkshift[1] .+ hval[h] .* bpkshift[2]) .+ 1)
                    p_instrument=pdf.(Binomial.(zbin-1,p_),z_dest-1)
                    xtran[:, x_dest+(z_dest-1)*xbin, h, j] .= p_mileage_transition .* p_instrument
                    lcdf1 = p_mileage_transition .+ lcdf1
                    lcdf2 = p_mileage_transition .* p_instrument .+ lcdf2
                    xtranc[:, x_dest+(z_dest-1)*xbin, h, j] .= xtranc[:, x_dest+(z_dest-1)*xbin, h, j] .+ lcdf2 # xtranc[:,x_dest] on the RHS here is completely superfluous
                end
            end
        end
    end
    

    xtran_1=zeros(precision, xbin_orig,phibin_dest,hbin,1)
    for h in 1:hbin
        xtran_1[:,:,h,1].=repeat(xtran[1,:,h,1]',xbin)
    end
    xtran=cat(xtran_1,xtran,dims=4)
    xtranc_1=cumsum(xtran[:,:,:,1],dims=2)
    xtranc=cat(xtranc_1,xtranc,dims=4)

    return xtran, xtranc
end

function xgrid_married(xtran) 
    xtran_married=zeros(precision, phibin_m_dest,phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1)
    xtranc_married_v1=copy(xtran_married)
    for x_m_orig in 1:xbin_orig
        for h_m in 1:hbin
            for x_f_orig in 1:xbin_orig
                for h_f in 1:hbin
                    for j_m in 1:J+1
                        for j_f in 1:J+1
                            for phi_m_dest in 1:phibin_dest
                                for phi_f_dest in 1:phibin_dest
                                    xtran_married[phi_m_dest,phi_f_dest,x_m_orig,x_f_orig,h_m,h_f,j_m,j_f]=xtran[x_m_orig,phi_m_dest,h_m,j_m]*xtran[x_f_orig,phi_f_dest,h_f,j_f]
                                end
                            end
                            # println("HERE")
                            xtranc_married_v1[:,:,x_m_orig,x_f_orig,h_m,h_f,j_m,j_f].=reshape(cumsum(xtran_married[:,:,x_m_orig,x_f_orig,h_m,h_f,j_m,j_f][:]),(phibin,phibin))
                        end
                    end
                end
            end
        end
    end

    xtranc_married_v2=cumsum(reshape(xtran_married,(phibin_m_dest*phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1)),dims=1) 
    xtranc_married_v2=reshape(xtranc_married_v2,(phibin_m_dest,phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1))
   
    return xtran_married, xtranc_married_v2
end

function xgrid_vect(b_fhat, bpkshift)
    # Turns out that this can be easily vectorized by taking at each step the difference  
    # (1-exp(b_fhat*(phival[i]-x))) - (1-exp(b_fhat*(phival[i-1]-x))). 
    # See my Goodnotes Misc Notes 3, search for ArcidiaconoMiller2011 transition process
    sign_1=zeros(precision, xbin,xbin)
    for x=1:xbin
        sign_1[x,x] = 1
        if x>1
            sign_1[x-1,x] = -1
        end
    end
    h_=reshape(hval,(1,1,hbin))
    xub = reshape(vcat(xval[2:end,1], precision.(Inf)),(1,:))
    b_fhat_=reshape(b_fhat,(1,1,1,J))
    
    # Calculate the differences once
    differences = xub .- xval
    # Logical matrix for the conditions
    condition_matrix = xub .> xval
    
    # Calculate the exponential term once
    exp_term = (1 .- exp.(-b_fhat_ .* h_ .*  differences)) .* condition_matrix
    
    # exp_term is now of size (phibin,phibin,hbin,J), and I need to matrix multiply each (:,:,h,j) by sign_1

    xtran_vect_1= zeros(precision, size(exp_term))
    for h in 1:hbin
        for j=1:J
            xtran_vect_1[:,:,h,j] =  exp_term[:,:,h,j] * sign_1
        end
    end

    xtran_1=zeros(precision, xbin,xbin,hbin,1)
    for h in 1:hbin
        xtran_1[:,:,h,1].=xtran_vect_1[1,:,h,1]'
    end
    xtran_vect_1=cat(xtran_1,xtran_vect_1,dims=4)
    # The following reshape is only used as a memo of the arrangement of the array
    xtran_vect_1=repeat(xtran_vect_1,1,zbin,1,1)

    # Probability distribution of the instrument
    p_= 1 ./(exp.(xval .* bpkshift[1] .+ reshape(hval,(1,1,:)) .* bpkshift[2]) .+ 1) 
    p_instrument=pdf.(Binomial.(zbin-1,p_), zstate' .-1) 
    p_instrument=reshape(p_instrument,(1,:,hbin,1))

    xtran_vect=p_instrument .* xtran_vect_1


    # Calculate the transition matrix of a couple. I assume that the transitions of each spouse is independent of the other
    xtran_vect_reshape=permutedims(xtran_vect,(2,1,3,4))
    xtran_married_vect=kron(xtran_vect_reshape[:],xtran_vect_reshape[:])
    xtran_married_vect=reshape(xtran_married_vect,(phibin_m_dest,xbin_m_orig,hbin_m,J_m+1,phibin_f_dest,xbin_f_orig,hbin_f,J_f+1)) # dimensions are (destination_male, origin_male, h_male, choice_male, destination_female, origin_female, h_female, choice_female)
    # rearrange dimensions as (destination_male, destination_female, origin_male, origin_female, h_male, h_female, choice_male, choice_female)
    xtran_married_vect=permutedims(xtran_married_vect,(1,5,2,6,3,7,4,8))

    

    # Calculate cumulative sums along columns
    xtranc_vect = cumsum(xtran_vect, dims=2)
    xtranc_married_vect=cumsum(reshape(xtran_married_vect,(phibin_m_dest*phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1)), dims=1) 
    xtranc_married_vect=reshape(xtranc_married_vect,(phibin_m_dest,phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1))

    return xtran_vect, xtranc_vect, xtran_married_vect, xtranc_married_vect
end
