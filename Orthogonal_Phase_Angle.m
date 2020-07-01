function [angPSF_pi2] = ...
            Orthogonal_Phase_Angle...
                (angPSF)
  

%% Calculate orthogonal angles to the PSF phase angles (+/-pi/2)___________
angPSF_pi2 = angPSF - pi/2;

%% Convert to complex representation_______________________________________
exp_iangPSF_pi2 = exp(1i*angPSF_pi2); 

%% Calculate phase of complex representation and its magnitude_____________
angPSF_pi2 = angle(exp_iangPSF_pi2);
mag_angPSF_pi2 = abs(angPSF_pi2);

%% Unwrap orthogonal angles________________________________________________
angPSF_pi2_p2pi = angPSF_pi2 + 2*pi;
            
angPSF_pi2(mag_angPSF_pi2>pi & angPSF_pi2<0) = ...
	angPSF_pi2_p2pi(mag_angPSF_pi2>pi & angPSF_pi2<0);

angPSF_pi2_m2pi = angPSF_pi2 - 2*pi;
            
angPSF_pi2(mag_angPSF_pi2>pi & angPSF_pi2>0) = ...
    angPSF_pi2_m2pi(mag_angPSF_pi2>pi & angPSF_pi2>0);

        
end