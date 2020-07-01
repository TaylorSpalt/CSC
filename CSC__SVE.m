function [arrayF_csc] = ...
            CSC__SVE...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb,SVE_sin_theta)


% All references to equations are from "Constrained Spectral Conditioning
%   for spatial sound level estimation", Spalt et al., 2016

%% Initialize______________________________________________________________     
%___Divide by 2/K in order to not have to scale other variables by 2/K
%   within the loop__________
aGamb_tK = aGamb/two_invK;   	
G_diag_tK = G_diag/two_invK;

%% Perform CSC, processing each channel individually_______________________
arrayF_csc = arrayF;            % Solution starts as initial dataset
for m = v_M
    %___Initialize reprocessing flag___
    reprocess = 0.5;

    %___Processing loop. The Fourier transform of the selected sensor is
    %   processed until no more reference sensors in the array are eligible
    %   to be used.________________________________________________________
    while reprocess > 0
        if reprocess == 0.5
            %___Start with unprocessed data___
            Gmx = G(:,m);	
            arrayF_m = arrayF(:,m);         	
        elseif reprocess == 1
            %___Use processed data once processed___
            Gmx = tK_cTarrayF*arrayF_m;             
        end

        %___Compute mag & phase of cross-spectra___
        aGmx = abs(Gmx);
        angGmx = angle(Gmx);

        %___Compute the cross-spectra cosine projections onto the 
        %   orthogonal phase vector______________
        cos_d_ang = cos(angPSF_pi2(:,m)-angGmx);

        %___Apply steering vector error constraints to cosine projection 
        %   term of Eq. (7)______________________
        [cos_d_ang_SVE] = ...
            SVE_Constraints...
                (cos_d_ang,SVE_sin_theta(:,m));
        
        %___Form weight vector magnitude via Eq. (7) using steering vector 
        %   error constraints (Eqs. (19-21) & Fig. 4)____
        Wb = aGmx.*cos_d_ang_SVE./G_diag;      

        %___Remove sensors which violate conditions of Eq. (9)___
        Wb(m) = NaN;            % Same mic cannot be reference                 
        Wb(abs(Wb)>1) = NaN;    % |Wb|>1 cannot be used     
        aWb = abs(Wb);

        %___Loop through ref sensors and process Fourier transform_____
        [aWb,arrayF_m,reprocess] = ...
            CSC_Core...
                (aWb,M,arrayF,Wb,exp_iangPSF_pi2,aGamb_tK,arrayF_m,...
                    reprocess,m);
        
        %___Once all ref sensors have been disqualified, terminate loop___
        if sum(isnan(aWb)) == M  
            %___As a safety check, only save result if its autospectral 
            %   mag is less than unprocessed autospectral mag__________
            if arrayF_m'*arrayF_m > G_diag_tK(m)  ||  reprocess == 0.5  
                arrayF_csc(:,m) = arrayF(:,m);
            else
                arrayF_csc(:,m) = arrayF_m;
            end
            
            %___Reset reprocessing flag___
            reprocess = 0;
        end
    end            
end

     
end