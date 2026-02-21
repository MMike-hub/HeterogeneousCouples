# Calculates FV terms given the data.
# Note: I have to calculate both continuation values for firt step and second step

function fvdata_RX1(CCP_,  xtran, xtran_married)
    global zeros_st1_m_part1, zeros_st1_m_part2, zeros_st1_f_part1, zeros_st1_f_part2, zeros_st2_cpl_m, zeros_st2_cpl_f, zeros_st2_m, zeros_st2_f
    

    ############################################################
    # Second-stage continuation values of married agents.
    ############################################################
    # Considering that in this model past spouses don't enter the agent's state and that there is a renewal action, I can build finite dependence in two periods by choosing to be single and replace engine in period t+1.
    # not renew --> single -->renew
    # renew --> single -->renew
    view1= @view CCP_.p0_st2_m[phistate_list[:,1],:,:,:]
    # The continuation value is E(-log p0_st1 - log p0_st2 +kappa(single|phi)) where the expectation is taken over the evolution of the observable state. See derivation on Lyx file. The kappa part cannot be precomputed since it's not nonparametrically identified. Instead, it will be a control variable in the structural regression.
    U_st2_married_m= -log.(CCP_.p0_st1_m_marriedb4) .- log.(reshape(view1,(phibin,hbin,NTypes,1,1,1,T2))) #The reshape here is to ensure that broadcasting subtraction works
    # Now take the expectation over the evolution of the couple's state.
    # remember that xtran_married is arranged as (phi_m_t,phi_f_t,x_m_tminus1,x_f_tminus1,h_m,h_f,a_m_tminus1,a_f_tminus1)
    # while U_st2_married_m is arranged (x_m_t,h_m,s_m,x_ftminus1_t,h_ftminus1,s_ftminus1,t) where ftminus1 indicates m's partner chosen at t minus 1 hence x_ftminus1_t indicates ftminus1's value of x at time t
    # So i need a permutedims and reshape to broadcast the multiplication correctly
    xtran_married_reshape=reshape(xtran_married,(phibin_m_dest, phibin_f_dest, xbin_m_orig, xbin_f_orig, hbin_m,hbin_f, J_m+1, J_f+1,1,1,1)) #add two singleton dimensions as placeholders for agents types, and one last dimension for time. 
    # Similarly, add placeholders dimensions to U_st2_married_m. Two of them, for male and female choices, plus two of them from origin state of the agent and his spouse (x_m_tminus1 and x_ftminus1_tminus1)
    U_st2_married_m_reshape=reshape(U_st2_married_m,(phibin_m, hbin_m, NTypes_m, phibin_sp_orig, hbin_sp_orig, NTypes_sp_orig, T2,1,1,1,1))
    # Then permute the dimensions so that they are in the same order: phi_m_t ,phi_ftminus1_t, x_m_tminus1, x_ftminus1_tminus1, h_m, h_ftminus1, s_m, s_ftminus1, a_m_tminus1, a_f_tminus1, t. 
    xtran_married_reshape=permutedims(xtran_married_reshape, (1,2,3,4,5,6,9,10,7,8,11))
    U_st2_married_m_reshape=permutedims(U_st2_married_m_reshape,(1,4,10,11,2,5,3,6,8,9,7))
    # Now multiply and sum across the destination state to compute the expectation
    EU_st2_married_m=sum(U_st2_married_m_reshape .* xtran_married_reshape,dims=(1,2))
    # Reshape to eliminate the two first singleton dimensions. Notice that because the spouse's type doesn't enter transitions, I enter NTypes_sp_orig instead of NTypes_f.
    EU_st2_married_m=reshape(EU_st2_married_m,(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_sp_orig,J_m+1,J_f+1,T2)) #(x_m_tminus1,x_ftminus1_tminus1,h_m,h_ftminus1,s_m,s_ftminus1,a_m_tminus1,a_f_tminus1,t) 
    # In the structural regression, I need to represent choice probabilities using the *difference* between the conditonal valuation functions between every choice and a baseline choice
    dEU_st2_married_m=Beta.*(EU_st2_married_m .-EU_st2_married_m[:,:,:,:,:,:,[1],[1],:])
    # Notice that EU_st2_married_m[:,:,:,:,:,:,[1],:,:]  and EU_st2_married_m[:,:,:,:,:,:,1,:,:] are different.The former will completely eliminate dimension 7, while the latter leaves dimension 7 as singleton.
    # Each value EU_st2_married_m[:,:,:,:,:,:,:,:,t] is the continuation value for an action chosen at t-1, so it will have to be summed to flow utilities of actions chosen at t-1 to form conditional valuation functions. For this reason, I need to decrease by one all values of t in EU_st2_married_m
    dEU_st2_married_m=cat(dEU_st2_married_m[:,:,:,:,:,:,:,:,2:end], zeros_st2_cpl_m,dims=9)
    # For each couple with characteristics x_m_tminus1,x_ftminus1_tminus1,h_m,h_ftminus1,s_m,s_ftminus1 making choice a_m_tminus1,a_f_tminus1 at time t-1, the part that is a function of CCPs of the expected value for the male in the couple from such choice is EU_st2_married_m.  
    # Now, rearrange as x_m,h_m,s_m,x_f,h_f,s_f,j_m,j_f,t, so that it remains arranged the same way as the couple CCPs
    dEU_st2_married_m=permutedims(dEU_st2_married_m,(1,3,5,2,4,6,7,8,9))


    #Now repeat to calculate the (part that is a function of CCPs of) continuation value of a married female
    # NOTICE THAT THE RESHAPE OF U_st2_married_f IS DIFFERENT THAN THAT OF U_st2_married_m
    view1= @view CCP_.p0_st2_f[phistate_list[:,1],:,:,:]
    U_st2_married_f= -log.(CCP_.p0_st1_f_marriedb4).- log.(reshape(view1,(phibin,hbin,NTypes,1,1,1,T2)))
    U_st2_married_f_reshape=reshape(U_st2_married_f,(phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2,1,1,1,1))
    U_st2_married_f_reshape=permutedims(U_st2_married_f_reshape,(4,1,11,10,5,2,6,3,9,8,7))
    EU_st2_married_f=sum(U_st2_married_f_reshape .* xtran_married_reshape,dims=(1,2))
    EU_st2_married_f=reshape(EU_st2_married_f,(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_sp_orig,NTypes_f,J_m+1,J_f+1,T2))
    dEU_st2_married_f=Beta.*(EU_st2_married_f .-EU_st2_married_f[:,:,:,:,:,:,[1],[1],:])
    dEU_st2_married_f=cat(dEU_st2_married_f[:,:,:,:,:,:,:,:,2:end], zeros_st2_cpl_f,dims=9)
    dEU_st2_married_f=permutedims(dEU_st2_married_f,(1,3,5,2,4,6,7,8,9))

    ##########################################################
    # Second-stage continuation values of single agents
    ##########################################################
    view1= @view CCP_.p0_st2_m[phistate_list[:,1],:,:,:]
    U_st2_single_m= -log.(CCP_.p0_st1_m_singleb4).-log.(view1)
    # xtran doesn't depend on unobserved type nor time, but U_st2_single_m does, so I add two placeholder dimensions to conform these array for element-wise multiplication.
    xtran_reshape=reshape(xtran,(xbin_m,phibin_m_dest,hbin_m,J_m+1,1,1)) #First dimension is origin x, second is destination phi, third is h, fourth is choice, fifth is unobs. type, sixth is time.
    # Now I permute the dimensions so that it's arranged in the order phi_t+1,x_t, h,s,choice at t,t
    xtran_reshape=permutedims(xtran_reshape,(2,1,3,5,4,6))
    #similarly, U_st2_single_m doesn't depend on the agent's choice, but xtran does. also xtran depends on previous-period x (x_origin) while U_st2_single_m doesn't. I add two placeholder dimensions
    U_st2_single_m_reshape=reshape(U_st2_single_m,(phibin_m,hbin_m,NTypes_m,T2,1,1))
    # Now I permute the dimensions so that it's arranged in the order x_t+1,x_t, h,s,choice at t,t
    U_st2_single_m_reshape=permutedims(U_st2_single_m_reshape,(1,6,2,3,5,4))
    EU_st2_single_m=sum(U_st2_single_m_reshape.*xtran_reshape,dims=1)
    #reshape to eliminate the first singleton dimension 
    # Now EU_st2_single_m is organized as (x_t,h,s,choice at t,t)
    EU_st2_single_m=reshape(EU_st2_single_m,(xbin_m,hbin_m, NTypes_m,J_m+1,T2))
    # Calculate the difference with respect to the baseline choice
    dEU_st2_single_m=Beta.*(EU_st2_single_m .- EU_st2_single_m[:,:,:,[1],:])
    # Since this is the CONTINUATION value for choices made at t, every choice at t should be matched at values of U_st2_single_m that use CCPs from t+1, that's why I have to shift all of dEU_st2_single_m down by one year. Obviously, the continuation values at T2 are zero.
    dEU_st2_single_m=cat(dEU_st2_single_m[:,:,:,:,2:end], zeros_st2_m,dims=5)

    view1= @view CCP_.p0_st2_f[phistate_list[:,1],:,:,:]
    U_st2_single_f=-log.(CCP_.p0_st1_f_singleb4).-log.(view1)
    U_st2_single_f_reshape=reshape(U_st2_single_f,(phibin_f,hbin_f,NTypes_f,T2,1,1))
    U_st2_single_f_reshape=permutedims(U_st2_single_f_reshape,(1,6,2,3,5,4))
    EU_st2_single_f=sum(U_st2_single_f_reshape.*xtran_reshape,dims=1)
    EU_st2_single_f=reshape(EU_st2_single_f,(xbin_f,hbin_f, NTypes_f,J_f+1,T2))
    dEU_st2_single_f=Beta.*(EU_st2_single_f .- EU_st2_single_f[:,:,:,[1],:])
    dEU_st2_single_f=cat(dEU_st2_single_f[:,:,:,:,2:end], zeros_st2_f,dims=5)

    ############################################################
    # First-stage continuation values
    ############################################################
    # In a model with transition costs in marriage, the only way to ensure finite dependence is to telescope conditional value functions until the next marital choice then choose to be single, then renew to reset mileage
    # marry a -------------------------> renew --> single --> renew
    # single (or marry any b \neq a) --> renew --> single --> renew
    view1= @view CCP_.p0_st2_m[phistate_list[:,1],:,:,:]
    U_st1_m_part1=-log.(CCP_.p0_st1_m_marriedb4).-log.(reshape(view1,(phibin_m,hbin_m,NTypes_m,1,1,1,T2)))
    # U_st1_m_part1 doesn't depend on the couple's second-stage choice, but xtran_married does. Also, U_st1_m_part1 doesn't depend on the origin x state of each spouse, but xtran_married does, so I need to add 4 placeholder dimensions for broadcasting multiplication. 
    U_st1_m_part1_reshape=reshape(U_st1_m_part1,(phibin_m, hbin_m, NTypes_m, phibin_sp_orig, hbin_sp_orig, NTypes_sp_orig, T2, 1,1,1,1))
    # Now I arrange the dimensions as phi_m_t1,phi_sp_t1,x_m_t, x_sp_t, h_m,h_sp,s_m,s_sp,a_m,a_sp,t
    U_st1_m_part1_reshape=permutedims(U_st1_m_part1_reshape,(1,4,8,9,2,5,3,6,10,11,7))
    # Similarly xtran_married doesn't depend on time nor the agents' unobservable type s and unobservable type s, so I need to add 3 placeholder dimensions to xtran_married
    xtran_married_reshape=reshape(xtran_married,(phibin_m_dest,phibin_f_dest,xbin_m_orig,xbin_f_orig,hbin_m,hbin_f,J_m+1,J_f+1,1,1,1))
    # now I rearrange the dimansions as phi_m_t1,phi_sp_t1, x_m_t, x_sp_t, h_m,h_sp,s_m,s_sp,a_m,a_sp,t
    xtran_married_reshape=permutedims(xtran_married_reshape,(1,2,3,4,5,6,9,10,7,8,11))
    EU_st1_m_part1=sum(U_st1_m_part1_reshape .* xtran_married_reshape, dims=(1,2))
    # Reshape to eliminate the first two singleton dimensions
    EU_st1_m_part1=reshape(EU_st1_m_part1,(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_m,NTypes_sp_orig,J_m+1,J_f+1,T2))
    # Shift all values by one year. Continuation values at t=T2 are zero
    EU_st1_m_part1=cat(EU_st1_m_part1[:,:,:,:,:,:,:,:,2:end], zeros_st1_m_part1,dims=9)
    # Rearrange dimensions as (x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t)
    EU_st1_m_part1=permutedims(EU_st1_m_part1,(1,3,5,2,4,6,7,8,9))
    
    view1= @view CCP_.p0_st2_m[phistate_list[:,1],:,:,:]
    U_st1_m_part2=-log.(CCP_.p0_st1_m_singleb4).-log.(view1)
    # U_st1_m_part2 doesn't depend on the agent's second-stage choice, but xtran does. Also, U_st1_m_part2 doesn't depend on the origin x state of the agent, but xtran does, so I need to add 2 placeholder dimensions for broadcasting multiplication
    U_st1_m_part2_reshape=reshape(U_st1_m_part2,(phibin_m,hbin_m,NTypes_m,T2,1,1))
    # Now, arrange the dimensions of U_st1_m_part2_reshape so that it's arranged as (x_m_dest,x_m_origin,h_m,s_m,a_m,t)
    U_st1_m_part2_reshape=permutedims(U_st1_m_part2_reshape,(1,5,2,3,6,4))
    # Multiply ONLY THE TRANSITION GIVEN THE BASELINE (RENEWAL) CHOICE and sum across the first dimension
    EU_st1_m_part2=sum(U_st1_m_part2_reshape.*xtran_reshape[:,:,:,:,[1],:],dims=1)
    # Reshape to eliminate the first and fifth stingleton dimensions
    EU_st1_m_part2=reshape(EU_st1_m_part2,(xbin_m,hbin_m,NTypes_m,T2))
    # Shift all values down by one year, continuation values at T2 are zero
    EU_st1_m_part2=cat(EU_st1_m_part2[:,:,:,2:end],zeros_st2_m[:,:,:,[1],:],dims=4)
    
    # The difference in continuation values is \sum_{a_m,a_f} p_st2_cpl .*[Beta .*(EU_st1_m_part1 .-EU_st1_m_part2) .-log.(p_rf_st2_cpl) .+log.(p0_st2_m)]
    # But to make all these terms comformable to element-wise subtraction, I need to add some placeholder dimensions and permute some dimensions. 
    # All of these terms need to be arranged as (x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t) 
    # EU_st1_m_part2 doesn't depend on x_f,h_f,s_f,a_m,a_f so I need to add 5 placeholder dimensions and permute
    EU_st1_m_part2=reshape(EU_st1_m_part2,(xbin_m,hbin_m,NTypes_m,T2,1,1,1,1,1))
    EU_st1_m_part2=permutedims(EU_st1_m_part2,(1,2,3,5,6,7,8,9,4))
    # p0_st2_cpl doesn't depend on a_m,a_f so I need to add 2 placeholder dimensions and permute
    p0_st2_cpl_reshape=reshape(CCP_.p0_st2_cpl,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2,1,1))
    p0_st2_cpl_reshape=permutedims(p0_st2_cpl_reshape,(1,2,3,4,5,6,8,9,7))
    # p0_st2_m does not depend on the female's characteristics nor any other actions a_m, a_f, so I add 3+2 placeholder dimensions
    p0_st2_m_reshape=reshape(CCP_.p0_st2_m,(xbin_m,hbin_m,NTypes_m,T2,1,1,1,1,1))
    p0_st2_m_reshape=permutedims(p0_st2_m_reshape,(1,2,3,5,6,7,8,9,4))
    dEU_st1_m=sum(CCP_.p_st2_cpl.*(Beta .*(EU_st1_m_part1 .-EU_st1_m_part2) .-log.(CCP_.p_st2_cpl) .+log.(p0_st2_m_reshape)),dims=(7,8))
    # dEU_st1_m=sum(CCP_.p_st2_cpl.*(Beta .*(EU_st1_m_part1 .-EU_st1_m_part2)  .+log.(p0_st2_m_reshape./CCP_.p_st2_cpl)),dims=(7,8))
    # Reshape to remove singleton dimensions 7 and 8
    dEU_st1_m=reshape(dEU_st1_m,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))


    # Now do the same for female. Keep in mind is that it changes the way I permute the dimensions of U_st1_f_part1_reshape and EU_st1_f_part2
    view1= @view CCP_.p0_st2_f[phistate_list[:,1],:,:,:]
    U_st1_f_part1=-log.(CCP_.p0_st1_f_marriedb4).-log.(reshape(view1,(phibin_f,hbin_f,NTypes_f,1,1,1,T2)))
    U_st1_f_part1_reshape=reshape(U_st1_f_part1,(phibin_f,hbin_f,NTypes_f,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2,1,1,1,1))
    U_st1_f_part1_reshape=permutedims(U_st1_f_part1_reshape,(4,1,9,8,5,2,6,3,11,10,7))
    EU_st1_f_part1=sum(U_st1_f_part1_reshape .* xtran_married_reshape, dims=(1,2))
    EU_st1_f_part1=reshape(EU_st1_f_part1,(xbin_m,xbin_f,hbin_m,hbin_f,NTypes_sp_orig,NTypes_f,J_m+1,J_f+1,T2))
    EU_st1_f_part1=cat(EU_st1_f_part1[:,:,:,:,:,:,:,:,2:end],zeros_st1_f_part1,dims=9)
    # Rearrange dimensions as (x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t)
    EU_st1_f_part1=permutedims(EU_st1_f_part1,(1,3,5,2,4,6,7,8,9))
    
    view1= @view CCP_.p0_st2_f[phistate_list[:,1],:,:,:]
    U_st1_f_part2=-log.(CCP_.p0_st1_f_singleb4).-log.(view1)
    U_st1_f_part2_reshape=reshape(U_st1_f_part2,(phibin_f,hbin_f,NTypes_f,T2,1,1))
    U_st1_f_part2_reshape=permutedims(U_st1_f_part2_reshape,(1,5,2,3,6,4))
    EU_st1_f_part2=sum(U_st1_f_part2_reshape.*xtran_reshape[:,:,:,:,[1],:],dims=1)
    EU_st1_f_part2=reshape(EU_st1_f_part2,(xbin_f,hbin_f,NTypes_f,T2))
    EU_st1_f_part2=cat(EU_st1_f_part2[:,:,:,2:end],zeros_st1_f_part2,dims=4)
    
    EU_st1_f_part2=reshape(EU_st1_f_part2,(xbin_f,hbin_f,NTypes_f,T2,1,1,1,1,1))
    EU_st1_f_part2=permutedims(EU_st1_f_part2,(5,6,7,1,2,3,8,9,4))
    p0_st2_cpl_reshape=reshape(CCP_.p0_st2_cpl,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2,1,1))
    p0_st2_cpl_reshape=permutedims(p0_st2_cpl_reshape,(1,2,3,4,5,6,8,9,7))
    p0_st2_f_reshape=reshape(CCP_.p0_st2_f,(xbin_f,hbin_f,NTypes_f,T2,1,1,1,1,1))
    p0_st2_f_reshape=permutedims(p0_st2_f_reshape,(5,6,7,1,2,3,8,9,4))
    dEU_st1_f=sum(CCP_.p_st2_cpl.*(Beta .*(EU_st1_f_part1 .-EU_st1_f_part2) .-log.(CCP_.p_st2_cpl) .+log.(p0_st2_f_reshape)),dims=(7,8))
    # dEU_st1_f=sum(CCP_.p_st2_cpl.*(Beta .*(EU_st1_f_part1 .-EU_st1_f_part2)  .+log.(p0_st2_f_reshape./CCP_.p_st2_cpl)),dims=(7,8))
    dEU_st1_f=reshape(dEU_st1_f,(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2))

    fvt1_RX1=FVT1_RX1(dEU_st1_f, dEU_st1_m, dEU_st2_single_f, dEU_st2_single_m, dEU_st2_married_f, dEU_st2_married_m)
    
    return fvt1_RX1
end

struct FVT1_RX1
    # This structure contains the continuation values for each point in the state space
    dEU_st1_f:: Union{Array, CuArray}
    dEU_st1_m:: Union{Array, CuArray}
    dEU_st2_single_f:: Union{Array, CuArray}
    dEU_st2_single_m:: Union{Array, CuArray}
    dEU_st2_married_f:: Union{Array, CuArray}
    dEU_st2_married_m:: Union{Array, CuArray}
end


