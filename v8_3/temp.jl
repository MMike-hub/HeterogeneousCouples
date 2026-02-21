 # Keep only T2 periods
        # Sampled_1=Sampled_1[sampled_1,1:T2]
        # Sampled_2=ifelse.(sampled_2.==0,0,Sampled_2[sampled_1,1:T2])
        State_1=State_1[sampled_1]
        State_2=ifelse.(sampled_2.==0,missing,State_2[sampled_1])
        State_sp_dest_1=State_sp_dest_1[sampled_1,1:T2]
        State_sp_dest_2=ifelse.(sampled_2.==0,missing,State_sp_dest_2[sampled_1,1:T2])
        FirstSampled_1= [findfirst(x -> x == 1, row) for row in eachrow(Sampled_1)][:]
        Y_1=ifelse.(Sampled_1.==0,missing,Y_1[sampled_1,1:T2])
        Y_sp_dest_1=ifelse.(Sampled_1.==0,missing,Y_sp_dest_1[sampled_1,1:T2])
        Y_2=ifelse.(Sampled_2.==0,missing,Y_2[sampled_1,1:T2])
        Y_sp_dest_2=ifelse.(Sampled_2.==0,missing,Y_sp_dest_2[sampled_1,1:T2])
        Phi_1=ifelse.(Sampled_1.==0,missing,Phi_1[sampled_1,1:T2,:])
        Phi_sp_dest_1=ifelse.(Sampled_1.==0,missing,Phi_sp_dest_1[sampled_1,1:T2,:])
        Phi_sp_orig_1=ifelse.(Sampled_1.==0,missing,Phi_sp_orig_1[sampled_1,1:T2,:])
        Phi_2=ifelse.(Sampled_2.==0,missing,Phi_2[sampled_1,1:T2,:])
        Phi_sp_dest_2=ifelse.(Sampled_2.==0,missing,Phi_sp_dest_2[sampled_1,1:T2,:])
        Phi_sp_orig_2=ifelse.(Sampled_2.==0,missing,Phi_sp_orig_2[sampled_1,1:T2,:])
        Phistate_1=ifelse.(Sampled_1.==0,missing,Phistate_1[sampled_1,1:T2])
        Phistate_sp_dest_1=ifelse.(Sampled_1.==0,missing,Phistate_sp_dest_1[sampled_1,1:T2])
        Phistate_sp_orig_1=ifelse.(Sampled_1.==0,missing,Phistate_sp_orig_1[sampled_1,1:T2])
        Phistate_2=ifelse.(Sampled_2.==0,missing,Phistate_2[sampled_1,1:T2])
        Phistate_sp_dest_2=ifelse.(Sampled_2.==0,missing,Phistate_sp_dest_2[sampled_1,1:T2])
        Phistate_sp_orig_2=ifelse.(Sampled_2.==0,missing,Phistate_sp_orig_2[sampled_1,1:T2])
        H_1=H_1[sampled_1]
        H_sp_dest_1=ifelse.(Sampled_1.==0,missing,H_sp_dest_1[sampled_1,1:T2])
        H_sp_orig_1=ifelse.(Sampled_1.==0,missing,H_sp_orig_1[sampled_1,1:T2])
        H_2=ifelse.(sampled_2.==0,missing,H_2[sampled_1])
        H_sp_dest_2=ifelse.(Sampled_2.==0,missing,H_sp_dest_2[sampled_1,1:T2])
        H_sp_orig_2=ifelse.(Sampled_2.==0,missing,H_sp_orig_2[sampled_1,1:T2])
        Hstate_1=Hstate_1[sampled_1]
        Hstate_sp_dest_1=ifelse.(Sampled_1.==0,missing,Hstate_sp_dest_1[sampled_1,1:T2])
        Hstate_sp_orig_1=ifelse.(Sampled_1.==0,missing,Hstate_sp_orig_1[sampled_1,1:T2])
        Hstate_2=ifelse.(sampled_2.==0,missing, Hstate_2[sampled_1])
        Hstate_sp_dest_2=ifelse.(Sampled_2.==0,missing,Hstate_sp_dest_2[sampled_1,1:T2])
        Hstate_sp_orig_2=ifelse.(Sampled_2.==0,missing,Hstate_sp_orig_2[sampled_1,1:T2])
        Sex_1=Sex_1[sampled_1]
        Sex_2=ifelse.(sampled_2.==0,missing,Sex_2[sampled_1])
        ID_1=ID_1[sampled_1]
        ID_2=ifelse.(sampled_2.==0,missing,ID_2[sampled_1])
        Married_1=ifelse.(Sampled_1.==0,missing,Married_1[sampled_1,1:T2])
        Married_2=ifelse.(Sampled_2.==0,missing,Married_2[sampled_1,1:T2])








        minimum(State_1_temp.===State_1_temp)
        minimum(State_2_temp.===State_2_temp)
        minimum(State_sp_dest_1_temp.===State_sp_dest_1_temp)
        minimum(State_sp_dest_2_temp.===State_sp_dest_2_temp)
        minimum(Y_1_temp.===Y_1_temp)
        minimum(Y_sp_dest_1_temp.===Y_sp_dest_1_temp)
        minimum(Y_2_temp.===Y_2_temp)
        minimum(Y_sp_dest_2_temp.===Y_sp_dest_2_temp)
        minimum(Phi_1_temp.===Phi_1_temp)
        minimum(Phi_sp_dest_1_temp.===Phi_sp_dest_1_temp)
        minimum(Phi_sp_orig_1_temp.===Phi_sp_orig_1_temp)
        minimum(Phi_2_temp.===Phi_2_temp)
        minimum(Phi_sp_dest_2_temp.===Phi_sp_dest_2_temp)
        minimum(Phi_sp_orig_2_temp.===Phi_sp_orig_2_temp)
        minimum(Phistate_1_temp.===Phistate_1_temp)
        minimum(Phistate_sp_dest_1_temp.===Phistate_sp_dest_1_temp)
        minimum(Phistate_sp_orig_1_temp.===Phistate_sp_orig_1_temp)
        minimum(Phistate_2_temp.===Phistate_2_temp)
        minimum(Phistate_sp_dest_2_temp.===Phistate_sp_dest_2_temp)
        minimum(Phistate_sp_orig_2_temp.===Phistate_sp_orig_2_temp)
        minimum(H_sp_dest_1_temp.===H_sp_dest_1_temp)
        minimum(H_sp_orig_1_temp.===H_sp_orig_1_temp)
        minimum(H_2_temp.===H_2_temp)
        minimum(H_sp_dest_2_temp.===H_sp_dest_2_temp)
        minimum(H_sp_orig_2_temp.===H_sp_orig_2_temp)
        minimum(Hstate_sp_dest_1_temp.===Hstate_sp_dest_1_temp)
        minimum(Hstate_sp_orig_1_temp.===Hstate_sp_orig_1_temp)
        minimum(Hstate_2_temp.===Hstate_2_temp)
        minimum(Hstate_sp_dest_2_temp.===Hstate_sp_dest_2_temp)
        minimum(Hstate_sp_orig_2_temp.===Hstate_sp_orig_2_temp)
        minimum(Sex_2_temp.===Sex_2_temp)
        minimum(ID_2_temp.===ID_2_temp)
        minimum(Married_1_temp.===Married_1_temp)
        minimum(Married_2_temp.===Married_2_temp)