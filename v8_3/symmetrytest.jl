test1=zeros(xbin,zbin,hbin,hbin,NTypes,NTypes,J+1,J+1,T2)
for x_m in 1:xbin
    for x_f in 1:xbin
        for h_m in 1:hbin
            for h_f in 1:hbin
                for s_m in 1:NTypes
                    for s_f in 1:NTypes
                        for j_m in 1:J+1
                            for j_f in 1:J+1
                                for t in 1:T2
                                    test1=abs.(p_st2_cpl[x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,t].==p_st2_cpl[x_f,x_m,h_f,h_m,s_f,s_m,j_f,j_m,t])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
minimum(test1)

3,3,1,1,2,1,3,1,7: test1= false, test2= false, test3=true, test4=true
3,3,2,2,2,1,2,3,7: test1= false, test2= false, test3=false, test4=false
3,3,2,2,2,1,3,1,1: test1= false, test2= false, test3=true, test4=true
3,3,2,2,2,2,1,2,6: test1= false, test2= false, test3=false, test4=false
3,3,2,2,2,1,3,3,6: test1= false, test2= false, test3=false, test4=false
3,3,2,2,2,1,1,1,3: test1= false, test2= false, test3=false, test4=false
util_married_fn(x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,sex,t)
x_m=3
x_f=3
h_m=1
h_f=1
s_m=2
s_f=1
j_m=3
j_f=1
t=7

x_m=3
x_f=3
h_m=1
h_f=1
s_m=1
s_f=2
j_m=1
j_f=3
t=7
vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,2,t]
sum(exp.(lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t].*vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,1,t] .+ (1 .-lambda[x_m,x_f,h_m,h_f,s_m,s_f,1,1,1,t]).*vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,:,:,2,t]))

test1=zeros(xbin,xbin,hbin,hbin,NTypes,NTypes,J+1,J+1,T2)
for x_m in 1:xbin
    for x_f in 1:xbin
        for h_m in 1:hbin
            for h_f in 1:hbin
                for s_m in 1:NTypes
                    for s_f in 1:NTypes
                        for j_m in 1:J+1
                            for j_f in 1:J+1
                                for t in 1:T2
                                    test1[x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,t]=vJ_married[x_m,x_f,h_m,h_f,s_m,s_f,j_m,j_f,1,t].==vJ_married[x_f,x_m,h_f,h_m,s_f,s_m,j_f,j_m,2,t]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
minimum(test1)

test1=zeros(xbin,xbin,hbin,hbin,NTypes,NTypes,T2)
test2=zeros(xbin,xbin,hbin,hbin,NTypes,NTypes,T2)
for x_m in 1:xbin
    for x_f in 1:xbin
        for h_m in 1:hbin
            for h_f in 1:hbin
                for s_m in 1:NTypes
                    for s_f in 1:NTypes
                        for t in 1:T2
                            test1[x_m,x_f,h_m,h_f,s_m,s_f,t]=abs.(FV_married[x_m,x_f,h_m,h_f,s_m,s_f,1,t].==FV_married[x_f,x_m,h_f,h_m,s_f,s_m,2,t])
                            test2[x_m,x_f,h_m,h_f,s_m,s_f,t]=abs.(FV_married[x_m,x_f,h_m,h_f,s_m,s_f,1,t].-FV_married[x_f,x_m,h_f,h_m,s_f,s_m,2,t])
                        end
                    end
                end
            end
        end
    end
end
minimum(test1)
maximum(test2)


