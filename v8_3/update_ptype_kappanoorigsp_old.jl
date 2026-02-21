# This file contains q update for the simplified case when the marriage transition costs do not depend on the origin spouse's characteristics

function comp(a,b) 
    """
    I'm just wrapping the identity function here to use "eachcol" operations... 
    """
    return a.===b 
end

function get_like(Like2, sp_ind)
    """
    This function returns the likelihood for a single time t for all indivudals.
    """
    type_p = Like2.*sp_ind
    replace!(type_p, 0.0 => 1.0)
    return prod(type_p, dims = 2)
end

function  update_ptype_old(Like, dt_CPU, Pi_cpl, Pi_single)
    """
    This function updates the probability that a couple is of a certain type
    Input 
        Like: vector of size 4N*T2*NTypes^2 that contains pi*p(a,phi)
        IS_sp: array that indicates the id of a person's partner in the sample
    Output
        PType: Probability that a Couple is of a given type (s_i, s_j). Size (4N,NTypes^2) 
    """
    global same_sp, First_t_with_sp, NotMarried

    if !RUN_ON_GPU
        # THIS METHOD DOES NOT WORK ON GPU
        if !@isdefined same_sp
            ID_sp_dest_=copy(dt.ID_sp_dest)
            same_sp = comp.(eachcol(ID_sp_dest_), Ref(ID_sp_dest_)) #This creates an array of T2 matrixes of size (2N*T2), entry (n,r) of matrix t indicates if the spouse at time t of guy n is the same as in time r.
        end
        Like2=copy(Array(Like))
        Like2 = reshape(Like2, 2*N, T2, NTypes^2) # Reshape the likelihood
        Like_sp_spells = get_like.(Ref(Like2), same_sp) #Get Like_sp_spells - this iterates across same_sp, the result is a vector of size T2 with 2N*1*NType^2 probs. (If the spouses were always the same, each of this T2 matrixes will be equal)
        Like_sp_spells=cat(Like_sp_spells...,dims=2) # Every element (n,t,s_cpl) contains the product of the likelihood contributions of all the periods that agent n spent married with their spouse at t, assuming the couple is of type (s,s_sp).
    else
        # THIS DOES
        if !@isdefined same_sp
            # No need to run this at every iteration
            ID_sp_dest_=copy(dt.ID_sp_dest)
            # ID_sp_dest_[dt.Married.===0].=Ref((0,0))
            same_sp = comp.(eachcol(ID_sp_dest_), Ref(ID_sp_dest_)) 
            same_sp=cat(same_sp..., dims=3)
            same_sp=CUDA.adapt(CuArray,same_sp)
        end
        Like2 = reshape(Like, 2*N, T2, 1, NTypes^2) 
        Like_sp_spells = Like2.* same_sp
        Like_sp_spells[Like_sp_spells.==0].=1
        Like_sp_spells=prod(Like_sp_spells,dims=2)
        Like_sp_spells=reshape(Like_sp_spells,(2*N,T2,NTypes^2))
    end
    # This reshape allows easy slice by agent type and spouse type
    # Along dimension 2 we have the two original spouses at the moment of sampling.
    Like_sp_spells=reshape(Like_sp_spells,(N,2,T2,NTypes,NTypes_sp))
    # Sanity check: verify that Like_sp_spells[i,2,t,:,:] is equal to one when agent i is single at t
    unique(Array(Like_sp_spells[:,2,:,:,:])[repeat(coalesce.(dt.Married[1:N,:],0),1,1,2,2).==0])

    # If n is single at time t, the element (n,t,s,s_sp) of Like_sp_spells contains the product of the likelihood contributions of all the periods that agent n spent single, assuming the agent is of type s, i.e. it should hold that (n,t,:,1)=(n,t,:,2) if the agent n is single at time t.
    # Sanity check.
    if !@isdefined NotMarried
        NotMarried= coalesce.(1 .-dt.Married,0)
        NotMarried=reshape(NotMarried,(N,2,T2))
        if RUN_ON_GPU
            NotMarried=CUDA.adapt(CuArray,NotMarried)
        end
    end
    check1= maximum(Like_sp_spells[:,:,:,:,1].*NotMarried .- Like_sp_spells[:,:,:,:,2].*NotMarried) > eps()
    if check1==true 
        println("ISSUE HERE: update_ptype_kappanoorigsp")
    end

    # Sanity check
    # FOR INDIVIDUALS i THAT WERE SAMPLED AS SINGLE, ALL VALUES OF Like_sp_spells_keep1[i,2,:,:,:] SHOULD BE EQUAL TO 1.
    # CHECK!
    selected_rows = findall(x -> x == 0, dt_CPU.MarriedHH[1:N])
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
    # VERIFY THAT q_s_sp__s==1/NTypes_sp WHEN THE AGENT IS SINGLE
    PType_sp_=q_s_sp__s[:]

    # q_sij, q_si
    # IMPORTANT: THIS FORMULA FOR q_s ONLY WORKS IF THE MARRIAGE TRANSITION COSTS DO NOT DEPEND ON THE ORIGIN SPOUSE
    if !@isdefined First_t_with_sp
        # First_t_with_sp[i,j] takes a value of 1 if the value of ID_sp_dest_[i,j] occurs in ID_sp_dest_[i,j] for the first time along the row ID_sp_dest_[i,:]
        First_t_with_sp= [[findfirst(Ref(ID_sp_dest_[i,j]) .=== ID_sp_dest_[i, :]).==j for j in 1:size(ID_sp_dest_,2)] for i in 1:size(ID_sp_dest_,1)]
        First_t_with_sp=vcat(transpose.(First_t_with_sp)...)
        First_t_with_sp=repeat(First_t_with_sp,1,1,NTypes,NTypes_sp)
        First_t_with_sp=reshape(First_t_with_sp,(N,2,T2,NTypes,NTypes_sp))
        if RUN_ON_GPU
            First_t_with_sp=CUDA.adapt(CuArray, First_t_with_sp)
        end
    end
    # Since Like_sp_spells contains at each (n,t,s,s_sp) the product of the likelihood contributions of ALL periods n spent with their spouse at t, I only need to take the first of such values for each spouse to construct the conditional likelihood needed for q_s
    Like_sp_spells_keep1=Like_sp_spells.*First_t_with_sp


    # Likelihood contribution of spells as single
    # Identify spells as single
    Like_single_spell=Like_sp_spells_keep1.*NotMarried
    # Only keep values when s_sp=1, since they are identical across all values of s_sp
    Like_single_spell=Like_single_spell[:,:,:,:,1]
    # Take the maximum values across time (There is at most one nonzero value along T2. They are all zero if the agent was never single)
    Like_single_spell=maximum(Like_single_spell,dims=3)
    # Replace zero values with ones
    # replace!(Like_single_spell, 0.0 => 1.0)
    Like_single_spell[Like_single_spell.==0].=1
    # Add a dimension for conformity to element-wise multiplication later and remove singleton time dimension
    Like_single_spell=reshape(Like_single_spell,(N,2,NTypes,1))

    # Likelihood contribution of spells married to original spouse
    Like_originalhh_spell=Like_sp_spells.*dt.OriginalSpouse_
    # Take the maximum values across time (There is at most one nonzero value along T2. They are all zero if the agent was single when sampled)
    Like_originalhh_spell=maximum(Like_originalhh_spell,dims=3)
    # Replace zero values with ones
    # replace!(Like_originalhh_spell, 0.0 => 1.0)
    Like_originalhh_spell[Like_originalhh_spell.==0].=1
    # Remove singleton time dimension
    Like_originalhh_spell=reshape(Like_originalhh_spell,(N,2,NTypes,NTypes_sp))

    # Now that I extracted from Like_sp_spells_keep1 the values of spells as single and as married to the original spouse, I can set those values to 0
    Like_sp_spells_keep1_v2=copy(Like_sp_spells_keep1)
    if !@isdefined NotMarried_v2
        NotMarried_v2=repeat(NotMarried,1,1,1,NTypes,NTypes_sp)
    end
    idx=findall(x -> x == 1, NotMarried_v2)
    Like_sp_spells_keep1_v2[idx].=0
    idx=findall(x -> x == 1, dt.OriginalSpouse_)
    Like_sp_spells_keep1_v2[idx].=0
    idx=findall(x -> x == 0, dt.insample)
    Like_sp_spells_keep1_v2[idx].=0

    Like_sp_spells_keep1_v2=sum(Like_sp_spells_keep1_v2,dims=5) #Integrate across the spouses' types
    # replace!(Like_sp_spells_keep1_v2, 0 => 1.0) # Replace 0 with ones since the next operation is a product
    Like_sp_spells_keep1_v2[Like_sp_spells_keep1_v2.==0].=1
    Like_newsp_spells=prod(Like_sp_spells_keep1_v2,dims=3) # Multiply likelihood contributions of all marriage spells across time
    # Remove singleton time dimension but leave singleton spouse type dimension
    Like_newsp_spells=reshape(Like_newsp_spells,(N,2,NTypes,1))

    # Multiply together the likelihoods of original spouses. Mind that you need to obtain all possible combination of types for the original couple
    # Along dimensions 3 and 4 of Like_newsp_spells_combined we have the types of the original couple
    Like_newsp_spells_combined=Like_newsp_spells[:,1,:,:] .* permutedims(Like_newsp_spells[:,2,:,:],(1,3,2))

    # Now combine the spells as single of both spouses in the original couple
    Like_single_spell_combined=Like_single_spell[:,1,:,:].*permutedims(Like_single_spell[:,2,:,:],(1,3,2))

    # Now combine the spells spent with the original spouse
    # I NEED TO DO THIS BECAUSE THE VALUES OF Like CORRESPONDING TO SUCH SPELLS ONLY CONTAIN THE FIRST-STAGE CHOICE PROBABILITY OF ONE OF THE TWO SPOUSES.
    # THEY BOTH CONTAIN THE SQUARE ROOT OF THE SECOND-STAGE CHOICE PROBABILITIES OF THE COUPLE, SO THAT THE PRODUCT WILL CONTAIN THE CORRECT LIKELIHOOD CONTRIBUTION
    Like_originalhh_spell_combined=Like_originalhh_spell[:,1,:,:] .* permutedims(Like_originalhh_spell[:,2,:,:],(1,3,2))
    
    
    Like_allspells=Like_newsp_spells_combined .* Like_single_spell_combined .* Like_originalhh_spell_combined


    # Apply Pi weights
    temp1=Like_allspells.*Pi_cpl
    temp2=Like_allspells.*Pi_single

    # Sanity check
    # Verify that temp2[:,:,1]==temp2[:,:,2] 
    minimum(temp2[:,:,1].==temp2[:,:,2] )

    q_sij=temp1./sum(temp1,dims=(2,3))
    
    q_si=temp2[:,:,1] ./ sum(temp2[:,:,1] , dims=2)
   
    #Sanity check
    if (maximum(sum(q_sij,dims=(2,3)))>1.0001) | (minimum(sum(q_sij,dims=(2,3)))<0) | (maximum( sum(q_si,dims=(2)))>1.0001) | (minimum(sum(q_si,dims=(2)))<0)
        println("issue in update_ptype: sum of q outside of [0,1]") 
    end
   
    q_sij_=reshape(q_sij,(N,1,NTypes,NTypes))
    q_sij_1=repeat(q_sij_,1,T2,1,1)
    q_sij_2=permutedims(q_sij_1,(1,2,4,3))
    q_sij_=vcat(q_sij_1,q_sij_2)
    
    q_si_=reshape(q_si,(N,1,NTypes))
    q_si_1=repeat(q_si_,1,T2,1,NTypes)
    q_si_= vcat(q_si_1,ones(size(q_si_1)))

    PType_=q_sij_[:].*dt.marriedhh .+ (1 .- dt.marriedhh).*q_si_[:]

    
    #Replace to missing when not in sample
    PType_ = ifelse.(dt.insample .== 1, PType_, 0.0)
    PType_sp_= ifelse.(dt.insample .== 1, PType_sp_, 0.0)
    

    PType=PType_.*PType_sp_

    return PType, q_sij, q_si
end

function  update_Pi_old(q_sij, q_si, dt)
    Pi_cpl_=sum(q_sij .* dt.MarriedHH[1:N],dims=1)./N
    Pi_single_=sum(q_si .* (1 .- dt.MarriedHH[1:N]),dims=1)./N
    return Pi_cpl_, Pi_single_
end