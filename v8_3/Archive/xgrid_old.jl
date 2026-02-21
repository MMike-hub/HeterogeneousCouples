function xgrid(theta, xval, j)
    n = length(xval)
    xub = vcat(xval[2:end], Inf)
    xtran = zeros(n, n)
    xtranc = zeros(n, n)
    lcdf = 0.0

    for i in 1:xbin
        # Consider the different discrete choices as partial engine replacements / maintenance.
        # The larger the value of j, the lower the probability that the next period mileage is high
        # i.e. the larger the j, the more thorough the maintenance / engine parts replacement
        xtran[:, i] .= (xub[i] .>= xval) .* (1 .- exp.(-theta .* j .* (xub[i] .- xval)) .- lcdf)
        lcdf = xtran[:, i] .+ lcdf
        xtranc[:, i] .= xtranc[:, i] .+ lcdf # xtranc[:,i] on the RHS here is completely superfluous
    end
    return xtran, xtranc
end

function xgrid_married(xtran)
    xtran_married=zeros(xbin_m_dest,xbin_f_dest,xbin_m_orig,xbin_f_orig,zbin_m,zbin_f,J_m+1,J_f+1)
    xtranc_married_v1=copy(xtran_married)
    for x_m_orig in 1:xbin
        for z_m in 1:zbin
            for x_f_orig in 1:xbin
                for z_f in 1:zbin
                    for j_m in 1:J+1
                        for j_f in 1:J+1
                            for x_m_dest in 1:xbin
                                for x_f_dest in 1:xbin
                                    xtran_married[x_m_dest,x_f_dest,x_m_orig,x_f_orig,z_m,z_f,j_m,j_f]=xtran[x_m_orig,x_m_dest,z_m,j_m]*xtran[x_f_orig,x_f_dest,z_f,j_f]
                                end
                            end
                            println("HERE")
                            xtranc_married_v1[:,:,x_m_orig,x_f_orig,z_m,z_f,j_m,j_f].=reshape(cumsum(xtran_married[:,:,x_m_orig,x_f_orig,z_m,z_f,j_m,j_f][:]),(xbin,xbin))
                        end
                    end
                end
            end
        end
    end
    xtranc_married_v2=cumsum(reshape(xtran_married,(xbin_m_dest*xbin_f_dest,xbin_m_orig,xbin_f_orig,zbin_m,zbin_f,J_m+1,J_f+1)),dims=1) 
    xtranc_married_v2=reshape(xtranc_married_v2,(xbin_m_dest,xbin_f_dest,xbin_m_orig,xbin_f_orig,zbin_m,zbin_f,J_m+1,J_f+1))
    return xtran_married, xtranc_married_v2
end

function xgrid_vect(lambda)
    # Turns out that this can be easily vectorized by taking at each step the difference  
    # (1-exp(lambda*(xval[i]-x))) - (1-exp(lambda*(xval[i-1]-x))). 
    # See my Goodnotes Misc Notes 3, search for ArcidiaconoMiller2011 transition process
    sign_1=zeros(xbin,xbin)
    for x=1:xbin
        sign_1[x,x] = 1
        if x>1
            sign_1[x-1,x] = -1
        end
    end
    z_=reshape(zval,(1,1,zbin))
    xub = reshape(vcat(xval[2:end], Inf),(1,:))
    lambda_=reshape(lambda,(1,1,1,J))
    
    # Calculate the differences once
    differences = xub .- xval
    # Logical matrix for the conditions
    condition_matrix = xub .> xval
    
    # Calculate the exponential term once
    exp_term = (1 .- exp.(-lambda_ .* z_ .*  differences)) .* condition_matrix
    
    # exp_term is now of size (xbin,xbin,zbin,J), and I need to matrix multiply each (:,:,z,j) by sign_1

    xtran_vect= zeros(size(exp_term))
    for z in 1:zbin
        for j=1:J
            xtran_vect[:,:,z,j] =  exp_term[:,:,z,j] * sign_1
        end
    end

    xtran_1=zeros(xbin,xbin,zbin,1)
    for z in 1:zbin
        xtran_1[:,:,z,1].=repeat(xtran_vect[1,:,z,1]',xbin)
    end
    xtran_vect=cat(xtran_1,xtran_vect,dims=4)
    # The following reshape is only used as a memo of the arrangement of the array
    xtran_vect=reshape(xtran_vect,(xbin_orig,xbin_dest,zbin,J+1)) 

    # Calculate the transition matrix of a couple. I assume that the transitions of each spouse is independent of the other
    xtran_vect_reshape=permutedims(xtran_vect,(2,1,3,4))
    xtran_married_vect=kron(xtran_vect_reshape[:],xtran_vect_reshape[:])
    xtran_married_vect=reshape(xtran_married_vect,(xbin_m_dest,xbin_m_orig,zbin_m,J_m+1,xbin_f_dest,xbin_f_orig,zbin_f,J_f+1)) # dimensions are (destination_male, origin_male, z_male, choice_male, destination_female, origin_female, z_female, choice_female)
    # rearrange dimensions as (destination_male, destination_female, origin_male, origin_female, z_male, z_female, choice_male, choice_female)
    xtran_married_vect=permutedims(xtran_married_vect,(1,5,2,6,3,7,4,8)) 
    

    # Calculate cumulative sums along columns
    xtranc_vect = cumsum(xtran_vect, dims=2)
    xtranc_vect=permutedims(xtranc_vect,(2,1,3,4))
    xtranc_vect=reshape(xtranc_vect,(xbin,xbin,zbin,J+1))
    xtranc_vect=permutedims(xtranc_vect,(2,1,3,4))
    xtranc_married_vect=cumsum(reshape(xtran_married_vect,(xbin_m_dest*xbin_f_dest,xbin_m_orig,xbin_f_orig,zbin_m,zbin_f,J_m+1,J_f+1)), dims=1) 
    xtranc_married_vect=reshape(xtranc_married_vect,(xbin_m_dest,xbin_f_dest,xbin_m_orig,xbin_f_orig,zbin_m,zbin_f,J_m+1,J_f+1))

    return xtran_vect, xtranc_vect, xtran_married_vect, xtranc_married_vect
end
