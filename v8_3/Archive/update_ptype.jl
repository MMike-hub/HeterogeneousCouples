function comp(a,b) 
    """
    I'm just wrapping the identity function here to use "eachcol" operations... 
    """
    return a.==b 
end

function get_like(Like2, sp_ind)
    """
    This function returns the likelihood for a single time t for all indivudals.
    """
    type_p = Like2.*sp_ind
    replace!(type_p, 0.0 => 1.0)
    return prod(type_p, dims = 2)
end
function  update_ptype(Like, ID_sp, insample)
    """
    This function updates the probability that a couple is of a certain type
    Input 
        Like: vector of size 4N*T*NTypes^2 that contains pi*p(a,phi)
        IS_sp: array that indicates the id of a person's partner in the sample
    Output
        PType: Probability that a Couple is of a given type (s_i, s_j). Size (4N,NTypes^2) 
    """
    same_sp = comp.(eachcol(ID_sp), Ref(ID_sp)) #This creates an array of T matrixes of size (2N*T), entry (n,r) of each matrix indicates if the spouse at time t of guy n is the same as in time r.
    Like2 = reshape(Like, 2*N, NTypes^2, T2) #Reshape the likelihood 
    replace!(Like2, missing => 1.0) #Replace the missing with 1.0 (multiplicative identity)
    
    PType = get_like.(Ref(permutedims(Like2,[1,3,2])), same_sp) #Get Ptype - this iterates across same_sp, the result is a vector of size T with 2N*NType^2 probs. (If the spouses were always the same, each of this T matrixes will be equal)
    PType = reshape(mapreduce(permutedims, hcat, reshape.(PType, 2*N*NTypes^2)), 2*N*T2*NTypes^2) #Change the dimensions to match all other vectors (e.g. Like)
    
    #Replace to missing when not in sample
    insample = reshape(repeat(dt.Insample, NTypes^2), 2*N*T2*NTypes^2)
    PType = ifelse.(insample .== 1, PType, missing)
    return PType
end
