function [aWb,arrayF_m,reprocess] = ...
            CSC_Core...
                (aWb,M,arrayF,Wb,exp_iangPSF_pi2,aGamb_tK,arrayF_m,...
                    reprocess,m)
 
%___Pre-subtraction Fourier transform magnitude___
mag_arrayF_m = arrayF_m'*arrayF_m;

%___Loop through ref sensors until Fourier transform is processed or all
%   ref sensors are disqualified_________________________________________
while sum(isnan(aWb)) < M
    %___Eq. (9)________________________________________
    [~,mp] = max(aWb);	% Ref is chosen as max|Wb|<=1

    %___Form weighted Fourier transform which is subtracted from
    %   Fourier transform of sensor being processed (Eq. (8))___
    arrayF_Wb = arrayF(:,mp)*Wb(mp)*exp_iangPSF_pi2(mp,m); 
    mag_arrayF_Wb = arrayF_Wb'*arrayF_Wb;

    %___If mag of weighted ref > ambient level (Eq. (13))...
    if mag_arrayF_Wb > aGamb_tK(m,mp)  
        % ...subtract weighted ref from input___________
        arrayF_m_p = arrayF_m - arrayF_Wb;     % Eq. (8)
        mag_arrayF_m_p = arrayF_m_p'*arrayF_m_p;

        % If the new magnitude < pre-subtraction magnitude (Eq. (12)), 
        %   save the updated Fourier transform___
        if mag_arrayF_m_p < mag_arrayF_m
            arrayF_m = arrayF_m_p;                 
            reprocess = 1;
            break;
        else
            %___Else, disqualify that ref sensor____
            aWb(mp) = NaN;
        end
    else
        %___Else, disqualify that ref sensor____
        aWb(mp) = NaN;
    end
end


end