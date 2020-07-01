function [gsxi,gsyi,gszi,Si] = ...
            Intragrid_Inputs...
                (s_intra,gsx,gsy,gsz)

            
%___Define grid spacing of intra-grid points in x-coord.___
if gsx < 1  
    gsxi = gsx/s_intra;          
    Sxi = s_intra + 1;
else
    gsxi = 1/s_intra;          
    Sxi = s_intra*gsx + 1;
end

%___Define grid spacing of intra-grid points in y-coord.___
if gsy < 1  
    gsyi = gsy/s_intra;          
    Syi = s_intra + 1;
else
    gsyi = 1/s_intra;         
    Syi = s_intra*gsy + 1;
end

%___Define grid spacing of intra-grid points in z-coord.___
if gsz < 1  
    gszi = gsz/s_intra;          
    Szi = s_intra + 1;
else
    gszi = 1/s_intra;          
    Szi = s_intra*gsz + 1;
end

%___Total # of intra-grid points___
Si = Sxi*Syi*Szi;   


end