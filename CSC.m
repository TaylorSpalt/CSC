function [Y_csc] = ...
            CSC...
                (e,cTe,G,G_diag,v_M,v_M_m1,M,DR,wY,v_S,z_S,v_G_diag,...
                    lamda,r_M_S,gsxi,gsyi,gszi,z_M_Si,dm,RE_switch,...
                    SVE_switch,RE_SVE_switch,arrayF,tK_cTarrayF,aGamb,...
                    two_invK,wG,sx,sy,sz,gsx,gsy,gsz,mx,my,mz,z_M_M,K,...
                    z_M,syn_data)


%% Use CSC to determine bform level at gmax found via MVDR or FDBF_________
Y_csc = z_S;
for s = v_S 
    %___Calculate inputs needed for CSC_________________________________
    [angPSF,angPSF_pi2,exp_iangPSF_pi2,SVE_sin_theta] = ...
        CSC_Inputs...
            (e(:,s),cTe(s,:),r_M_S(:,s),sx(s),sy(s),sz(s),gsxi,gsyi,...
                gszi,z_M_Si,dm,SVE_switch,RE_SVE_switch,v_M_m1,lamda,...
                gsx,gsy,gsz,mx,my,mz,M,v_G_diag,z_M_M,syn_data,v_M);

    %___Based on constraint selection, perform CSC_________________________
    if RE_switch == 1
        %___Random error constraints______________________________________
        [arrayF_csc] = ...
            CSC__RE...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb,K,angPSF,z_M);
    elseif SVE_switch == 1
        %___Steering vector error constraints_____________________________
        [arrayF_csc] = ...
            CSC__SVE...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb,SVE_sin_theta);
    elseif RE_SVE_switch == 1
        %___Random & steering vector error constraints____________________
        [arrayF_csc] = ...
            CSC__RE_SVE...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb,K,SVE_sin_theta,angPSF,z_M);
    else
        %___No constraints________________________________________________
        [arrayF_csc] = ...
            CSC__Unconstrained...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb);
    end
    
    %__Form CSM from CSC output_________________________
    G_csc = two_invK*ctranspose(arrayF_csc)*arrayF_csc;
    
    %___Modify CSM diagonal accordingly_____________
    if DR == 0
        G_csc(v_G_diag) = max(0,real(diag(G_csc)));
    end

    %__Beamform to targeted grid point with CSC-CSM_________
    Y_csc(s) = max(0,real(cTe(s,:)*(wG.*G_csc)*e(:,s)))*wY; 
end


end