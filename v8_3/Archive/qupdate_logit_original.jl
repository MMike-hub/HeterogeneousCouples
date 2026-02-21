function qupdate_logit_original(j,binit,binits,base, intcondX,lp)

if j > 1
    intcond_objective = b -> intcond(b,base, intcondX)
    if j<50
        fun=TwiceDifferentiable(intcond_objective,binits; autodiff = :forward)
        result = optimize(fun, binits, BFGS(), o2)
    end
    if j>=50
        fun=TwiceDifferentiable(intcond_objective,binit; autodiff = :forward)
        result = optimize(fun, binit, BFGS(), o2)
    end
    binit = Optim.minimizer(result)
end

templ = intcond(binit, base, intcondX)
lp = vcat(lp, templ)

PType = intcondP(binit, base, intcondX)
PType = repeat(PType, T2)
PType = reshape(PType, (N * T2 * 2))

    return PType, binit, lp
end