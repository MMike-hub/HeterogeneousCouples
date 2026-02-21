function intcond(b, like, X)
    global jjj

    jjj = b
    U1 = X * b
    p = exp.(U1) ./ (1 .+ exp.(U1))
    p = hcat(p, 1 .- p)

    Like = -sum(log.(max.(sum(p .* like, dims=2), eps()))) / 200

    return Like
end
