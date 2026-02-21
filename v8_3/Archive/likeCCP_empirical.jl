#  This is essentially the same as wlogit, in the sense that the logit 
#  probabilities are calculated over the space of polynomials of x and z,
#  i.e. the matrix of regressors is xx. But it's the same as likeCCP in that
#  instead of calculating the overall
#  "empirical" log likelihood, it calculates the vector of likelihood
#  contributions of each observed choice. The first half of the Like output
#  vector is going to be the CCPs of the observed choices conditional on
#  s=0, the second half is going to be the CCPs of the same observed
#  choices, in the exact same order, but conditional on s=1

function likeCCP_empirical(b1, y2, xx)
    dim=size(xx)[1]
    dim_xx=size(xx)[2]
    b1=reshape(b1,(dim_xx,J))
    b1=hcat(zeros(dim_xx,1),b1)
    U1 = xx * b1

    idx = CartesianIndex.(collect(1:dim), y2[:].+1)

    Like = exp.(U1[idx]) ./ sum(exp.(U1),dims=2)
    return Like
end
