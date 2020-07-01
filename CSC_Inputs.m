function [angPSF,angPSF_pi2,exp_iangPSF_pi2,SVE_sin_theta] = ...
            CSC_Inputs...
                (e_s,cTe_s,r_M_s,sx_s,sy_s,sz_s,gsxi,gsyi,gszi,z_M_Si,...
                    dm,SVE_switch,SVE_RE_switch,v_M_m1,lamda,gsx,gsy,...
                    gsz,mx,my,mz,M,v_G_diag,z_M_M,syn_data,v_M)

    
%% Calculate Point Spread Function from each grid point to all sensors_____
PSF = e_s*cTe_s;

%___Phase angle of the PSF___
angPSF = angle(PSF); 

%% Calculate orthogonal angles to the PSF phase angles_____________________
[angPSF_pi2] = ...
    Orthogonal_Phase_Angle...
        (angPSF); 

%___Convert orthogonal angles to complex representation____
exp_iangPSF_pi2 = exp(1i*angPSF_pi2);

%% Calculate the uncertainty in steering vector phase angle between sensor 
%   pairs due to steering vector errors____________________________________       
if SVE_switch == 1  ||  SVE_RE_switch == 1        
    [SVE_sin_theta] = ...
        SVE_Phase_Uncertainty...
            (sx_s,sy_s,sz_s,gsx,gsy,gsz,mx,my,mz,M,v_G_diag,lamda,...
                z_M_Si,gsxi,gsyi,gszi,dm,z_M_M,r_M_s,v_M_m1,syn_data,v_M);
else
    SVE_sin_theta = NaN;
end


end