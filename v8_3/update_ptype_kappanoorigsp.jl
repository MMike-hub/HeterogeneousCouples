# This file contains q update for the simplified case when the marriage transition costs do not depend on the origin spouse's characteristics



function get_like(Like2, sp_ind)
    """
    This function returns the likelihood for a single time t for all indivudals.
    """
    type_p = Like2.*sp_ind
    replace!(type_p, 0.0 => 1.0)
    return prod(type_p, dims = 2)
end

function  update_ptype(Like, dt, Pi_cpl, Pi_single_m, Pi_single_f)
    """
    This function updates the probability that a couple is of a certain type
    Input 
        Like: vector of size 4N*T2*NTypes^2 that contains pi*p(a,phi)
        IS_sp: array that indicates the id of a person's partner in the sample
    Output
        PType: Probability that a Couple is of a given type (s_i, s_j). Size (4N,NTypes^2) 
    """
    
    # The get_like function was heavy on memory and pretty slow.
    # I need to use Float64 for this operation because the likelihoods calculated here are very small and often default to 0 if I use lower levels of precision
    Like2 = Float64.(reshape(Like, 2*N, T2, 1, NTypes^2) )
    Like_sp_spells = Like2.* dt.same_sp
    Like_sp_spells[Like_sp_spells.==0].=1
    Like_sp_spells=prod(Like_sp_spells,dims=2)
    Like_sp_spells=reshape(Like_sp_spells,(2*N,T2,NTypes^2))
    
    # This reshape and rearrangement allows easy slice by gender, type, and spouse type
    # Along dimension 2 we have the two household members.
    Like_sp_spells=reshape(Like_sp_spells,(N,2,T2,NTypes,NTypes_sp))
    # Sanity check: verify that Like_sp_spells[i,2,:,:,:] is equal to one when agent i was single at sampling
    unique(Array(Like_sp_spells[:,2,:,:,:])[Array(repeat(dt.MarriedHH[1:N],1,T2,2,2).==0)])
    # Sanity check: If n is single at time t, then Like_sp_spells[n,t,s,s_sp] contains the product of the likelihood contributions of all the periods that agent n spent single, assuming the agent is of type s, i.e. it should hold that (n,t,:,1)=(n,t,:,2) if the agent n is single at time t.
    # Sanity check.
    check1= maximum(Like_sp_spells[:,:,:,:,1].*dt.NotMarried .- Like_sp_spells[:,:,:,:,2].*dt.NotMarried) > eps()
    if check1==true 
        println("ISSUE HERE: update_ptype_kappanoorigsp")
    end
    # Sanity check
    # FOR INDIVIDUALS i THAT WERE SAMPLED AS SINGLE, ALL VALUES OF Like_sp_spells[i,2,:,:,:] SHOULD BE EQUAL TO 1.
    # CHECK!
    selected_rows = findall(x -> x == 0, dt.MarriedHH[1:N])
    if isempty(selected_rows)==false
        # Slice the Like_sp_spells_keep1 array using selected_rows and get the minimum
        result = minimum(Like_sp_spells[selected_rows, 2, :, :, :])
        if result<1
            println("Issue here in update_ptype: all these values should be 1")
        end
    end

    # q_{s_sp|s}
    # IMPORTANT: THIS FORMULA FOR q_{s_sp|s} ONLY WORKS IF THE MARRIAGE TRANSITION COSTS DO NOT DEPEND ON THE ORIGIN SPOUSE
    q_s_sp__s=Like_sp_spells ./ sum(Like_sp_spells,dims=5)
    # VERIFY THAT q_s_sp__s==1 WHEN THE AGENT IS SINGLE
    PType_sp_=q_s_sp__s[:]

    # q_sij, q_si
    # IMPORTANT: THIS FORMULA FOR q_s ONLY WORKS IF THE MARRIAGE TRANSITION COSTS DO NOT DEPEND ON THE ORIGIN SPOUSE
    # Since Like_sp_spells contains at each (n,t,s,s_sp) the product of the likelihood contributions of ALL periods n spent with their spouse at t, I only need to take the first of such values for each spouse to construct the conditional likelihood needed for q_s
    Like_sp_spells_keep1=Like_sp_spells.*dt.First_t_with_sp

    # Now, sort Like_sp_spells so that Like_sp_spells[:,1,:,:] contains males and Like_sp_spells[:,2,:,:] contains females
    Like_sp_spells_keep1m=Like_sp_spells_keep1[:,[1],:,:,:].*(dt.Sex_[1:N,1].==1) .+ Like_sp_spells_keep1[:,[2],:,:,:].*(dt.Sex_[1:N,1].==2)
    Like_sp_spells_keep1f=Like_sp_spells_keep1[:,[1],:,:,:].*(dt.Sex_[1:N,1].==2) .+ Like_sp_spells_keep1[:,[2],:,:,:].*(dt.Sex_[1:N,1].==1)
    Like_sp_spells_keep1_sort=cat(Like_sp_spells_keep1m, Like_sp_spells_keep1f, dims=2) 
 
    # Likelihood contribution of spells as single
    # Identify spells as single
    Like_single_spell=Like_sp_spells_keep1_sort.*dt.NotMarried_sort
    # Only keep values when s_sp=1, since they are identical across all values of s_sp
    # Sanity check
    # minimum(Like_single_spell[:,:,:,:,1].==Like_single_spell[:,:,:,:,2])
    Like_single_spell=Like_single_spell[:,:,:,:,[1]]
    # Take the maximum values across time (There is at most one nonzero value along T2. They are all zero if the agent was never single)
    Like_single_spell=maximum(Like_single_spell,dims=3)
    # Replace zero values with ones
    # replace!(Like_single_spell, 0.0 => 1.0)
    Like_single_spell[Like_single_spell.==0].=1
    # Remove singleton time dimension
    Like_single_spell=reshape(Like_single_spell,(N,2,NTypes,1))

    # Likelihood contribution of spells married to original spouse
    Like_originalhh_spell=Like_sp_spells_keep1_sort .* dt.OriginalSpouse_sort
    # Take the maximum values across time (There is at most one nonzero value along T2. They are all zero if the agent was single when sampled)
    Like_originalhh_spell=maximum(Like_originalhh_spell,dims=3)
    # Replace zero values with ones
    # replace!(Like_originalhh_spell, 0.0 => 1.0)
    Like_originalhh_spell[Like_originalhh_spell.==0].=1
    # Remove singleton time dimension
    Like_originalhh_spell=reshape(Like_originalhh_spell,(N,NSex,NTypes,NTypes_sp))

    # Now that I extracted from Like_sp_spells_keep1_sort the values of spells as single and as married to the original spouse, I can set those values to 1
   
    Like_sp_spells_keep1_sort[dt.idx_NotMarried_sort].=0
    Like_sp_spells_keep1_sort[dt.idx_OriginalSpouse_sort].=0
    Like_sp_spells_keep1_sort[dt.idx_insample_sort].=0

    Like_sp_spells_keep1_sort=sum(Like_sp_spells_keep1_sort,dims=5) #Integrate across the spouses' types
    # replace!(Like_sp_spells_keep1_sort, 0 => 1.0) # Replace 0 with ones since the next operation is a product
    Like_sp_spells_keep1_sort[Like_sp_spells_keep1_sort.==0].=1
    Like_newsp_spells=prod(Like_sp_spells_keep1_sort,dims=3) # Multiply likelihood contributions of all marriage spells across time
    # if true
    #     # when working with low precision levels, prod(prod(Like_sp_spells_keep1_sort,dims=3)) can return zero values even if all elements of Like_sp_spells_keep1_sort are greater than zero.
    #     Like_newsp_spells[Like_newsp_spells.==0].=nextfloat(zero(precision))
    # end
    # Remove singleton time dimension but leave singleton spouse type dimension
    Like_newsp_spells=reshape(Like_newsp_spells,(N,2,NTypes,1))

    # Multiply together the likelihoods of original spouses. Mind that you need to obtain all possible combination of types for the original couple
    # Along dimensions 3 and 4 of Like_newsp_spells_combined we have the types of the original couple
    Like_newsp_spells_combined=Like_newsp_spells[:,1,:,:] .* permutedims(Like_newsp_spells[:,2,:,:],(1,3,2))
    # if true
    #     Like_newsp_spells_combined[Like_newsp_spells_combined.==0].=nextfloat(zero(precision))
    # end

    # Now combine the spells as single of both spouses in the original couple
    Like_single_spell_combined=Like_single_spell[:,1,:,:].*permutedims(Like_single_spell[:,2,:,:],(1,3,2))

    # Now combine the spells spent with the original spouse
    # I NEED TO DO THIS BECAUSE THE VALUES OF Like CORRESPONDING TO SUCH SPELLS ONLY CONTAIN THE FIRST-STAGE CHOICE PROBABILITY OF ONE OF THE TWO SPOUSES.
    # THEY BOTH CONTAIN THE SQUARE ROOT OF THE SECOND-STAGE CHOICE PROBABILITIES OF THE COUPLE, SO THAT THE PRODUCT WILL CONTAIN THE CORRECT LIKELIHOOD CONTRIBUTION
    Like_originalhh_spell_combined=Like_originalhh_spell[:,1,:,:] .* permutedims(Like_originalhh_spell[:,2,:,:],(1,3,2))
    
    # along the first dimension, we have originally sampled households, along the second we have male type, along the third we have female type
    Like_allspells=Like_newsp_spells_combined .* Like_single_spell_combined .* Like_originalhh_spell_combined

    # Apply Pi weights
    temp1=Like_allspells.*Pi_cpl
    temp2_m=Like_allspells.* Pi_single_m
    temp2_f=Like_allspells.* reshape(Pi_single_f,(1,1,NTypes))

    # Sanity check
    # Verify that temp2[:,:,1]==temp2[:,:,2] for single males and temp2[:,1,:]==temp2[:,2,:] for single females
    minimum(Like_allspells[(dt.MarriedHH[1:N].==0).*(dt.Sex_[1:N,1].==1),:,1].==Like_allspells[(dt.MarriedHH[1:N].==0).*(dt.Sex_[1:N,1].==1),:,2] )
    minimum(Like_allspells[(dt.MarriedHH[1:N].==0).*(dt.Sex_[1:N,1].==2),1,:].==Like_allspells[(dt.MarriedHH[1:N].==0).*(dt.Sex_[1:N,1].==2),2,:] )
    unique(Array(Like_allspells[(dt.Sex_[1:N,1].==1).*(dt.MarriedHH[1:N].==0),:,:]))

    q_sij=temp1./sum(temp1,dims=(2,3))
    
    q_si_m=(temp2_m[:,:,1] ./ sum(temp2_m[:,:,1] , dims=2)).*(dt.Sex_[1:N,1].==1)
    q_si_f=(temp2_f[:,1,:] ./ sum(temp2_f[:,1,:] , dims=2)).*(dt.Sex_[1:N,1].==2)
    q_si=q_si_m .+ q_si_f

    if maximum(dt.marriedhh)==false
        q_sij[:,:,:].=0
    end
    if minimum(dt.marriedhh)==true
        q_si[:,:].=0
    end

    #Sanity check
    if (maximum(sum(q_sij,dims=(2,3)))>1.0001) | (minimum(sum(q_sij,dims=(2,3)))<0) | (maximum( sum(q_si,dims=(2)))>1.0001) | (minimum(sum(q_si,dims=(2)))<0)
        println("issue in update_ptype: sum of q outside of [0,1]") 
    end
    
    # Right now q_sij is organized as (household, male_type, female_type) and I need to convert this back into the (individual, individual_type, spouse_type) format 
    q_sij_m=reshape(q_sij,(N,1,NTypes_m,NTypes_f))
    q_sij_m=repeat(q_sij_m,1,T2,1,1)
    q_sij_f=permutedims(q_sij_m,(1,2,4,3))
    q_sij_1=q_sij_m.*(dt.Sex_[1:N].==1).+q_sij_f.*(dt.Sex_[1:N].==2)
    q_sij_2=permutedims(q_sij_1,(1,2,4,3))
    q_sij_=vcat(q_sij_1,q_sij_2)
    
    q_si_=reshape(q_si,(N,1,NTypes))
    q_si_1=repeat(q_si_,1,T2,1,NTypes)
    q_si_= vcat(q_si_1,ones(precision, size(q_si_1)))

    PType_=q_sij_[:].*dt.marriedhh .+ (1 .- dt.marriedhh) .* q_si_[:]

    # Replace to missing when not in sample
    # PType_ = ifelse.(dt.insample .== 1, PType_, 0.0)
    # PType_sp_= ifelse.(dt.insample .== 1, PType_sp_, 0.0)
    PType_=PType_.*dt.insample
    PType_sp_=PType_sp_.*dt.insample

    PType=PType_.*PType_sp_

    return PType, q_sij, q_si
end

function  update_Pi(q_sij, q_si, dt)
    Pi_cpl_=sum(q_sij .* dt.MarriedHH[1:N],dims=1)./N
    Pi_single_m=sum(q_si .* (1 .- dt.MarriedHH[1:N]) .* (dt.Sex_[1:N].==1),dims=1)./N
    Pi_single_f=sum(q_si .* (1 .- dt.MarriedHH[1:N]) .* (dt.Sex_[1:N].==2),dims=1)./N
    return Pi_cpl_, Pi_single_m, Pi_single_f
end