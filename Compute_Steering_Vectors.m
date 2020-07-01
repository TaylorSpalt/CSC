function [e,cTe,r_M_S] = ...
            Compute_Steering_Vectors...
                (sx,sy,sz,mx,my,mz,ik,z_M_S,v_M,syn_data)
      
     
% Formulation II from E. Sarradj, "Three-dimensional acoustic source 
%   mapping with different beamforming steering vector formulations", 
%   Adv. Acoust. Vib. (2012).

%% Calculate Euclidean distances between sensors and grid points___________
r_M_S = z_M_S; 
for m = v_M
    r_M_S(m,:) = (sx(:)-mx(m)).^2 + (sy(:)-my(m)).^2 + (sz(:)-mz(m)).^2;
end
r_M_S = sqrt(r_M_S);

%% Determine center array sensor___________________________________________
[~,mc] = min(sqrt((mx-0).^2+(my-0).^2));    

%% Calculate raw steering vectors__________________________________________
if syn_data == 1
    e0 = exp(-ik*r_M_S);    % Flip wavenumber sign if synthetic data   
else
    e0 = exp(ik*r_M_S);
end

%% Normalize raw steering vectors by ratio of distance between sensor and
%   center array sensor____________________________________________________
e = z_M_S; 
for m = v_M
    e(m,:) = (r_M_S(m,:)./r_M_S(mc,:)).*e0(m,:);
end
cTe = ctranspose(e);


end