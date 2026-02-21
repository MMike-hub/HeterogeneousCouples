# This functions generates an array that contains the values of x x_sp_orig, z, z_sp_orig, s, s_sp_orig, and time dummies along the first dimension. The other dimensions index all points in the space of x x_sp_orig, z, z_sp_orig, s, s_sp_orig, and time dummies.

# THIS IS VERY SIMILAR TO RX1 GENERATED IN reducedform_state_original BUT THE VALUES IN DIMENSION 1 ARE NOW THE REGRESSORS FOR THE STRUCTURAL LOGIT, NOT THE REDUCED FORM ONE
function f_Xstruct_grid()
    # First stage.
    dim_Xstruct_st1_marriedb4_1=Int(length(bk)) # Number of regressors in the structural regression
    Xstruct_st1_marriedb4=zeros(precision, dim_Xstruct_st1_marriedb4_1,phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,phibin_sp_orig,hbin_sp_orig,NTypes_sp_orig,T2)
    for phi in 1:phibin
        for h in 1:hbin
            for s in 1:NTypes
                for phi_sp_orig in 1:phibin_sp_orig
                    for h_sp_orig in 1:hbin_sp_orig
                        for s_sp_orig in 1:NTypes_sp_orig
                            for phi_sp_dest in 1:phibin_sp_dest
                                for h_sp_dest in 1:hbin_sp_dest
                                    for s_sp_dest in 1:NTypes_sp_dest
                                        for t=1:T2
                                            vect=[]
                                            # vect=vcat(vect,phival[phi_sp_dest,1])
                                            # vect=vcat(vect,phival[phi_sp_dest,2])
                                            # vect=vcat(vect, hval[h_sp_dest])
                                            # vect=vcat(vect, s_sp_dest)
                                            # vect=vcat(vect, phival[phi,1])
                                            # vect=vcat(vect, phival[phi,2])
                                            # vect=vcat(vect, hval[h])
                                            # vect=vcat(vect, s)
                                            # vect=vcat(vect, 1) # Dummy for previously married
                                            # vect=vcat(vect, phival[phi_sp_orig,1])
                                            # vect=vcat(vect, phival[phi_sp_orig,2])
                                            # vect=vcat(vect, hval[h_sp_orig])
                                            # vect=vcat(vect, s_sp_orig)
                                            vect=vcat(vect,1)
                                            # tdummy=zeros(T2-1)
                                            # if t>1
                                            #     tdummy[t-1]=1
                                            # end
                                            # vect=vcat(vect, tdummy)
                                            Xstruct_st1_marriedb4[:,phi_sp_dest,h_sp_dest,s_sp_dest,phi,h,s,phi_sp_orig,h_sp_orig,s_sp_orig,t].=vect

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    dim_Xstruct_st1_singleb4_1=Int(length(bk))
    Xstruct_st1_singleb4=zeros(precision, dim_Xstruct_st1_singleb4_1,phibin_sp_dest,hbin_sp_dest,NTypes_sp_dest,phibin,hbin,NTypes,T2)
    for phi in 1:phibin
        for h in 1:hbin
            for s in 1:NTypes
                for phi_sp_dest in 1:phibin_sp_dest
                    for h_sp_dest in 1:hbin_sp_dest
                        for s_sp_dest in 1:NTypes_sp_dest
                            for t=1:T2
                                vect=[]
                        
                                # vect=vcat(vect,phival[phi_sp_dest,1])
                                # vect=vcat(vect,phival[phi_sp_dest,2])
                                # vect=vcat(vect, hval[h_sp_dest])
                                # vect=vcat(vect, s_sp_dest)
                                # vect=vcat(vect, phival[phi,1])
                                # vect=vcat(vect, phival[phi,2])
                                # vect=vcat(vect, hval[h])
                                # vect=vcat(vect, s)
                                # The following are placeholders, in case I want to multiply the same vector of parameters by Xstruct_st1_singleb4 and Xstruct_st1_marriedb4
                                # vect=vcat(vect, 0) # dummy for previously married
                                # vect=vcat(vect, 0) # Placeholder for phival[phi_sp_orig,1]
                                # vect=vcat(vect, 0) # Placeholder for phival[phi_sp_orig,2]
                                # vect=vcat(vect, 0) # Placeholder for hval[h_sp_orig]
                                # vect=vcat(vect, 0) # Placeholder for s_sp_orig
                                vect=vcat(vect,1)
                                # tdummy=zeros(T2-1)
                                # if t>1
                                #     tdummy[t-1]=1
                                # end
                                # vect=vcat(vect, tdummy)

                                Xstruct_st1_singleb4[:,phi_sp_dest,h_sp_dest,s_sp_dest,phi,h,s,t].=vect
                            end
                        end
                    end
                end
            end
        end
    end

    # Second stage.
    # dim_Xstruct_st2_single_1=Int(4+(T2-1))
    dim_Xstruct_st2_single_1=Int(4)
    Xstruct_st2_single=zeros(precision, dim_Xstruct_st2_single_1,xbin,hbin,NTypes,J+1,T2)
    for x in 1:xbin
        for h in 1:hbin
            for s in 1:NTypes
                for t=1:T2
                    for a in 1:J+1
                        # Establish normalization of flow utilities
                        if a ==1
                            Xstruct_st2_single[:,x,h,s,a,t].=0
                        else
                            vect=[]
                            vect=vcat(vect, 1)
                            vect=vcat(vect, xval[x])
                            # vect=vcat(vect, hval[h])
                            vect=vcat(vect, (s==2))
                            vect=vcat(vect, a)

                            # tdummy=zeros(T2-1)
                            # if t>1
                            #     tdummy[t-1]=1
                            # end
                            # vect=vcat(vect, tdummy)

                            Xstruct_st2_single[:,x,h,s,a,t].=vect
                        end
                    end
                end
            end
        end
    end

    # The dimensions of Xstruct_st2_m_married after the first one index x, x_sp, z, z_sp, s, s_sp, choice of agent, choice of spouse and the remaining are time dummies
    # dim_Xstruct_st2_m_married_1=Int(4+(T2-1))
    dim_Xstruct_st2_m_married_1=Int(4)
    Xstruct_st2_m_married=zeros(precision, dim_Xstruct_st2_m_married_1,xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2)
    for x_m in 1:xbin
        for h_m in 1:hbin
            for s_m in 1:NTypes
                for x_f in 1:xbin
                    for h_f in 1:hbin
                        for s_f in 1:NTypes
                            for a_m in 1:J_m+1
                                for a_f in 1:J_f+1
                                    for t=1:T2
                                        # Establish normalization of flow utilities
                                        if (a_m==1) 
                                            Xstruct_st2_m_married[:,x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t].=0
                                        else
                                            lambda=1
                                            # lambda=Lambda[x_m,h_m,s_m,x_f,h_f,s_f,1,1,t]
                                            vect=[]
                                            vect=vcat(vect, lambda * 1)
                                            vect=vcat(vect, lambda * xval[x_m])
                                            # vect=vcat(vect, lambda * hval[h_m])
                                            vect=vcat(vect, lambda * (s_m==2))
                                            vect=vcat(vect, lambda * a_m)
                                            # vect=vcat(vect, lambda * xval[x_f])
                                            # vect=vcat(vect, lambda * hval[h_f])
                                            # vect=vcat(vect, lambda * s_f)
                                            # vect=vcat(vect, lambda * a_f)
                                            # tdummy=zeros(T2-1)
                                            # if t>1
                                            #     tdummy[t-1]=1
                                            # end
                                            # vect=vcat(vect, tdummy)

                                            Xstruct_st2_m_married[:,x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t].=vect
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # The dimensions of Xstruct_st2_m_married after the first one index x, x_sp, z, z_sp, s, s_sp, choice of agent, choice of spouse and the remaining are time dummies
    # dim_Xstruct_st2_f_married_1=Int(4+(T2-1))
    dim_Xstruct_st2_f_married_1=Int(4)
    Xstruct_st2_f_married=zeros(precision, dim_Xstruct_st2_f_married_1,xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,J_m+1,J_f+1,T2)
    for x_m in 1:xbin
        for h_m in 1:hbin
            for s_m in 1:NTypes
                for x_f in 1:xbin
                    for h_f in 1:hbin
                        for s_f in 1:NTypes
                            for a_m in 1:J_m+1
                                for a_f in 1:J_f+1
                                    for t=1:T2
                                        # Establish normalization of flow utilities
                                        if (a_f==1)
                                            Xstruct_st2_f_married[:,x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t].=0
                                        else
                                            lambda=1
                                            # lambda=1-Lambda[x_m,h_m,s_m,x_f,h_f,s_f,1,1,t]
                                            vect=[]
                                            vect=vcat(vect, lambda * 1)
                                            vect=vcat(vect, lambda * xval[x_f])
                                            # vect=vcat(vect, lambda * hval[h_m])
                                            vect=vcat(vect, lambda * (s_f==2))
                                            vect=vcat(vect, lambda * a_f)
                                            # vect=vcat(vect, lambda * xval[x_f])
                                            # vect=vcat(vect, lambda * hval[h_f])
                                            # vect=vcat(vect, lambda * s_f)
                                            # vect=vcat(vect, lambda * a_f)
                                            # tdummy=zeros(T2-1)
                                            # if t>1
                                            #     tdummy[t-1]=1
                                            # end
                                            # vect=vcat(vect, tdummy)

                                            Xstruct_st2_f_married[:,x_m,h_m,s_m,x_f,h_f,s_f,a_m,a_f,t].=vect
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Linear indices to be used in CCP_structural. Convert from stage 2 space to stage 1 space
    dim_state_st2=(xbin_m,hbin_m,NTypes_m,xbin_f,hbin_f,NTypes_f,T2)
    idx_state_st2=reshape(collect(1:prod(dim_state_st2)), dim_state_st2)
    idx_st2_into_st1_f=repeat(idx_state_st2,zbin,1,1,zbin,1,1,1) #THIS IS ZBIN NOT HBIN
    idx_st2_into_st1_m=permutedims(idx_st2_into_st1_f,(4,5,6,1,2,3,7))

    Xstruct_grid=XstructGrid(Xstruct_st1_marriedb4, Xstruct_st1_singleb4, Xstruct_st2_m_married, Xstruct_st2_f_married, Xstruct_st2_single, idx_st2_into_st1_f, idx_st2_into_st1_m,  dim_Xstruct_st1_marriedb4_1, dim_Xstruct_st1_singleb4_1, dim_Xstruct_st2_single_1, dim_Xstruct_st2_m_married_1, dim_Xstruct_st2_f_married_1)

    return Xstruct_grid
end


struct XstructGrid
    Xstruct_st1_marriedb4:: Array
    Xstruct_st1_singleb4:: Array
    Xstruct_st2_m_married:: Array
    Xstruct_st2_f_married:: Array
    Xstruct_st2_single:: Array
    idx_st2_into_st1_f:: Array
    idx_st2_into_st1_m:: Array

    dim_Xstruct_st1_marriedb4_1:: Int
    dim_Xstruct_st1_singleb4_1:: Int
    dim_Xstruct_st2_single_1:: Int
    dim_Xstruct_st2_m_married_1:: Int
    dim_Xstruct_st2_f_married_1:: Int
end

# This function splits the structural parameter vector bccp_.
# This needs to be adapted to (1) the structure of the arrays in Xstruct_grid above and (2) restrictions such as assuming that some utility parameters are the same for men and women
function bccp_split(bccp_, Xstruct_grid)
    dim_Xstruct_st1_marriedb4_1=Xstruct_grid.dim_Xstruct_st1_marriedb4_1
    dim_Xstruct_st1_singleb4_1=Xstruct_grid.dim_Xstruct_st1_singleb4_1
    dim_Xstruct_st2_single_1= Xstruct_grid.dim_Xstruct_st2_single_1
    dim_Xstruct_st2_m_married_1=Xstruct_grid.dim_Xstruct_st2_m_married_1
    dim_Xstruct_st2_f_married_1=Xstruct_grid.dim_Xstruct_st2_f_married_1

    bccp_st1_m_marriedb4 = get_slice(bccp_, 1, dim_Xstruct_st1_marriedb4_1)
    bccp_st1_m_singleb4 = bccp_st1_m_marriedb4
    bccp_st1_f_marriedb4 = bccp_st1_m_marriedb4
    bccp_st1_f_singleb4 = bccp_st1_m_marriedb4
    bccp_st2_m_single = get_slice(bccp_, (dim_Xstruct_st1_marriedb4_1)+1, dim_Xstruct_st2_single_1)
    bccp_st2_f_single = bccp_st2_m_single
    bccp_st2_m_married = bccp_st2_m_single
    bccp_st2_f_married = bccp_st2_m_single
    return bccp_st1_m_marriedb4, bccp_st1_m_singleb4, bccp_st1_f_marriedb4, bccp_st1_f_singleb4, bccp_st2_m_single, bccp_st2_f_single, bccp_st2_m_married, bccp_st2_f_married
end