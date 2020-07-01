function [SVE_sin_theta] = ...
            SVE_Phase_Uncertainty...
                (sx_s,sy_s,sz_s,gsx,gsy,gsz,mx,my,mz,M,v_G_diag,lamda,...
                    z_M_Si,gsxi,gsyi,gszi,dm,z_M_M,r_M_s,v_M_m1,...
                    syn_data,v_M)
                        

%% Define intragrid point scan volume around center of grid point being 
%   targeted_______________________________________________________________
x_start = -gsx/2 + sx_s;
x_end = gsx/2 + sx_s;
y_start = -gsy/2 + sy_s;
y_end = gsy/2 + sy_s;
z_start = -gsz/2 + sz_s;
z_end = gsz/2 + sz_s;
[sxi,syi,szi] = ...
    meshgrid...
        (x_start:gsxi:x_end,y_start:gsyi:y_end,z_start:gszi:z_end);
        
%% Find max additional delta |distance| measured between the sensors and the
%   intragrid points, when compared to the nominal distance between the 
%   sensors and the grid point being targeted______________________________  

%___Initialize variables____
dmax_M_si = z_M_M;

%___If data is computer-generated, no uncertainty in sensor positions
%   exists______________________________________________________________
if syn_data == 1
    %___Distance from each intragrid point to each sensor____
    r_M_Si = z_M_Si;
    for m = v_M
        r_M_Si(m,:) = (sxi(:)-mx(m)).^2 + ...
                        (syi(:)-my(m)).^2 + ...
                    	(szi(:)-mz(m)).^2;
    end   
    r_M_Si = sqrt(r_M_Si);

    %___Calculate increase in delta |distance| due to intragrid points___
    for m = v_M_m1
        for mp = m+1:M
            %___|Delta| distance to grid point targeted___
            d_r_s = abs(r_M_s(m)-r_M_s(mp));
            
            %___|Delta| distance to intragrid points___
            d_r_si = abs(r_M_Si(m,:)-r_M_Si(mp,:));
            
            %___Maximum additional delta |distance| from all of the 
            %   intragrid points_________________
            dmax_M_si(m,mp) = max(d_r_si-d_r_s);
        end
    end
    
else %-------------------------------------------------------------------
    
    %___Real world sensor positions will always have some error_________
    for nmxp = [-dm,dm]
        mxp = mx + nmxp;
        for nmyp = [-dm,dm]
            myp = my + nmyp;
            for nmzp = [-dm,dm]
                mzp = mz + nmzp;
                
                %___Distance from each intragrid point to each sensor___
                r_M_Si = z_M_Si;
                for m = v_M
                    r_M_Si(m,:) = (sxi(:)-mxp(m)).^2 + ...
                                    (syi(:)-myp(m)).^2 + ...
                                    (szi(:)-mzp(m)).^2;
                end   
                r_M_Si = sqrt(r_M_Si);  

                %___Calculate increase in delta |distance| due to 
                %   intragrid points_____________________________________
                dmax_Mp_si = z_M_M;
                for m = v_M_m1
                    for mp = m+1:M
                        %___Nominal |difference| in distance___
                        d_r_s = abs(r_M_s(m)-r_M_s(mp));

                        %___Intragrid |difference| in distance___
                        d_r_si = abs(r_M_Si(m,:)-r_M_Si(mp,:));

                        %___Maximum additional delta |distance| from all 
                        %   of the intragrid points___________
                        dmax_Mp_si(m,mp) = max(d_r_si-d_r_s);
                    end
                end

                %__Update with max values______________
                dmax_M_si = max(dmax_Mp_si,dmax_M_si);
            end
        end
    end
end

%___Fill in lower half of matrix______________
dmax_M_si = dmax_M_si + transpose(dmax_M_si);

%___Set any decreases in delta |distance| to 0_____
dmax_M_si = max(0,dmax_M_si);

%___Convert delta distances into a % of wavelength_____
dmax_M_si = dmax_M_si/lamda;

%___Distances greater than 1/4 wavelength (i.e., pi/2) cannot be used 
%   because beyond pi/2 ambiguity exists in the phase (Eq. (19))______
dmax_M_si(dmax_M_si>=0.25) = NaN;

%___Auto pairs cannot be used_____
dmax_M_si(v_G_diag) = NaN;

%___Convert to radians and take sin in order to project source components
%	onto orthogonal vector used for cancellation (Fig. 4)______
SVE_sin_theta = sin(dmax_M_si*2*pi);


end