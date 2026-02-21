

function gmm_struct_closure(Y_dis, X_dis,Z_dis, PType_dis,W, i=1)
    # Remember that the orthogonality condition is x_j'(d_jnt-p_j(x)) for all j, j' in {1...J}, so I have to create one orth condition for each combination of state-variable and choice residual
    dim=size(X_dis)[1]
    dim_X=size(X_dis)[2]
    # X_v1=reshape(X_dis,(dim,dim_X*J))
    X_v1=[]
    Z_v1=reshape(repeat(Z_dis,1,J),(dim,dim_X*J*J))
    X_v2=reshape(permutedims(X_dis,(2,1,3)),(dim_X,dim*J))'
    Y_indicator=Int.(zeros((dim,J)))
    for y =1:J
        Y_indicator[:,y]=(Y_dis.==y)
    end
    Y_indicator=repeat(kron(Y_indicator,ones((1,dim_X))),1,J)
    # Use this if optimizing with Optim
    return bccp -> gmm_struct(bccp,Y_indicator,X_v1,X_v2,Z_v1, PType_dis,W, dim,dim_X, i)[i]
    # Use this if optimizing with Optimization
    # return (bccp,_) -> gmm_struct(bccp,Y_indicator,X_v1,X_v2,Z_v1, PType_dis,W, idx, dim,dim_X, )[i]
end


function gmm_struct(bccp,Y_indicator,X_v1,X_v2,Z_v1, PType_dis,W,  dim,dim_X, i)
    U1 = X_v2 * bccp
    U1=reshape(U1,(dim,J))
    U1=hcat(zeros(dim,1), U1)
    U1_=kron(U1[:,2:end],ones(1,dim_X))
    U1_=kron(ones(1,J),U1_)


    res=Y_indicator .- (exp.(U1_) ./ sum(exp.(U1), dims=2))
    prod_ =Z_v1 .* res 

    prod_weighted=PType_dis .* prod_

    # Sum across types
    # dim_s=Int(dim/NTypes)
    # prod_weighted=reshape(prod_weighted,(dim_s,NTypes,dim_X*J*J))
    # prod_weighted=permutedims(prod_weighted,(1,3,2))
    # prod_weighted=sum(prod_weighted,dims=3)[:,:,1]

    dim_s=dim
    if i in [1,  3]
        moment=sum(prod_weighted,dims=1)[:]./ dim_s
        if i==1
            obj=sum(moment' * W * moment)
        else
            obj=[]
        end
    else
        moment=[]
        obj=[]
    end

    if i ==2
        Σ=prod_weighted' * prod_weighted ./ dim_s
    else
        Σ=[]
    end    

    return obj, Σ, moment, prod_weighted
end