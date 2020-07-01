function [arrayF_csc] = ...
            CSC__RE...
                (M,v_M,G,arrayF,tK_cTarrayF,two_invK,G_diag,angPSF_pi2,...
                    exp_iangPSF_pi2,aGamb,K,angPSF,z_M)
                

% All references to equations are from "Constrained Spectral Conditioning
%   for spatial sound level estimation", Spalt et al., 2016

%% Initialize______________________________________________________________     
%___Divide by 2/K in order to not have to scale other variables by 2/K
%   within the loop__________
aGamb_tK = aGamb/two_invK;   	
G_diag_tK = G_diag/two_invK;

%___Variables needed for random error constraints___
TK = 2*K;
pi2 = pi/2;
sqrt2 = sqrt(2);

%___Calculate then apply normalized error in autospectral magnitude 
%   (Eqs. (14,16))___________
E_Gxx = 1/sqrt(K);
G_diag_RE = G_diag*(1+E_Gxx);

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
        
        %___Apply random error constraints to cross-spectral magnitude and
        %   cosine projection terms of Eq. (7)___________________________
        [aGmx_RE,cos_d_ang_RE] =...
            RE_Constraints...
                (aGmx,G_diag,K,TK,pi2,angPSF_pi2(:,m),angGmx,sqrt2,m,...
                    angPSF(:,m),z_M,cos_d_ang);
        
        %___Form weight vector magnitude via Eq. (7) using random error 
        %   constraints (Eqs. (14-17) & Fig. 2)____
        Wb = aGmx_RE.*cos_d_ang_RE./G_diag_RE;   
        
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