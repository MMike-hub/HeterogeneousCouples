# This calculates the likelihood of the sample using the input matrix of CCPs P0.

function likeCCP_empirical_nonpar(y2,Xstate,Zstate, P0)
    T2_=size(Xstate)[2]
    z2 = (Zstate .- 1) .* xbin .+ 1
    XZstate = Int.(z2 .+ Xstate .- 1)
    XZstate=reshape(XZstate,(N*T2_))
    XZstate=repeat(XZstate,NTypes)
    Tstate=repeat(reshape(collect(1:T2_),(1,:)),N)
    Tstate=reshape(Tstate,(N*T2_))
    Tstate=repeat(Tstate,NTypes)
    Sstate=Int.(vcat(ones(N*T2_), 2 .* ones(N*T2_)))
    P0_empiric=zeros(size(XZstate))
    for i=1:size(XZstate)[1]
        P0_empiric[i]=P0[XZstate[i],Sstate[i],Tstate[i]]
    end
    
    Like = y2 .* P0_empiric .+  (1 .-y2).*(1 .- P0_empiric)
    return Like
end
