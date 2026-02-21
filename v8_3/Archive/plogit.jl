function plogit(b, X)
    dim_X=size(X)[2]
    b=reshape(b,(dim_X,J))
    b=hcat(zeros(dim_X,1),b)
    U1 = X * b

    Like = 1 ./ sum(exp.(U1),dims=2)
    return Like
end
