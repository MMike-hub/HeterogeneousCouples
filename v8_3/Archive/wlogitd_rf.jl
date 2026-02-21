function wlogitd_rf(b1_, xx, PType, y2_indicator, xx_)
    dim_xx=size(xx)[2]
    b1_=reshape(b1_,(dim_xx,J))
    b1_=hcat(zeros(dim_xx,1),b1_)
    U1 = xx * b1_

    # compute the gradient
    # We have as many derivative as parameters, and we have dim_xx*J parameters. dg is of size (N*T2*Ntype,dim_xx*J)
    # The first dim_xx columns of dg represent the derivatives with respect to the parameters corresponding to choice j=1,
    # the following dim_xx columns are the derivatives w.r.t. the choice j=2 and so on.
    U1_=kron(U1[:,2:end],ones(1,dim_xx))
    dg = PType .* ( y2_indicator .- (exp.(U1_) ./ sum(exp.(U1),dims=2))).* xx_
    return dg
end