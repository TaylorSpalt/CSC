function [Y_true] = ...
            True_Source_Distribution...
                (S,nsource,cTe,e,Gs,wY,xsource,ysource,zsource,sx,sy,...
                    sz,wG,z_M_M,z_S,v_S,v_nsource)
        
        
%% Reshape scanning grid___________________________________________________
sx_re = reshape(sx,S,1);
sy_re = reshape(sy,S,1);
sz_re = reshape(sz,S,1);

%% Determine closest grid point to each source_____________________________
r_S_source = zeros(S,nsource);
min_r_S_source = zeros(nsource,1);
for source = v_nsource
    r_S_source(:,source) = sqrt((xsource(source)-sx_re).^2 + ...
                             	(ysource(source)-sy_re).^2 + ...
                            	(zsource(source)-sz_re).^2);
    [~,min_r_S_source(source)] = min(r_S_source(:,source));
end

%% Sum source CSMs which fall closest to each grid point then beamform
%   with the summed CSMs at that grid point to obtain the true source 
%   value there____________________________________________________________
Y_true = z_S;
for s = v_S
    %___For each grid point, determine which source CSM(s) fall closest
    %   to it, then sum all that do_________
    G_sum = z_M_M;
    for source = v_nsource
        if min_r_S_source(source) == s     
            G_sum = G_sum + Gs(:,:,source);  
        end
    end
    
    %____Compute the frequency-domain beamform using the summed CSM_____
    Y_true(s) = cTe(s,:)*(wG.*G_sum)*e(:,s);   
end

%___Apply FDBF normalization constant and set any grid points without
%   sources to NaN________________
Y_true = max(0,real(Y_true))*wY;
Y_true(Y_true==0) = NaN;


end