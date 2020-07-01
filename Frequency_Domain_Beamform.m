function [Y_fdbf] = ...
            Frequency_Domain_Beamform...
                (e,cTe,G,wY,v_S,z_S,wG)
                

%___Loop through grid points and calculate FDBF___
G = wG.*G;
Y_fdbf = z_S;
for g = v_S
    Y_fdbf(g) = cTe(g,:)*G*e(:,g);
end

%___Discard negative values as they don't have physical meaning, then
%   normalize beamform by normalization constant____
Y_fdbf = max(0,real(Y_fdbf))*wY;


end