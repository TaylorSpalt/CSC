function [G,G_diag,Gs,nsource,xsource,ysource,zsource,v_nsource,...
            arrayF,tK_cTarrayF,aGamb] = ...
                Generate_FFT_Array_Data...
                    (K,SNR,spert,sz0,c,two_invK,nFFT,ftrue,mx,my,mz,M,...
                        v_G_diag,v_M,v_K)

        
%% Specify point source location(s), inches________________________________
[xsource,ysource,zsource,nsource] = ...
    Source_Coordinates...
        (spert,sz0);
v_nsource = 1:nsource;

%% Define relative source pressure magnitudes______________________________
P_mag0 = [1,2,4,8]';

%% Calculate geometric variables___________________________________________
%___Distance from source(s) to sensors_______
r_M_nsource = zeros(M,nsource);
for source = v_nsource
    r_M_nsource(:,source) = ...
        abs(xsource(source)-mx(:)).^2 + ...
        abs(ysource(source)-my(:)).^2 + ...
        abs(zsource(source)-mz(:)).^2;  
end
r_M_nsource = sqrt(r_M_nsource);

%___Calculate average distance to sensors for each source___
r_M_nsource_avg = mean(r_M_nsource,1)';

%___Adjust source pressure magnitudes to account for distance between 
%   source(s) and sensors so distribution is even across the sensors_____
P_mag0 = P_mag0.*(r_M_nsource_avg/max(r_M_nsource_avg));

%% Create source pressure data_____________________________________________
P = zeros(K,nsource);
for source = v_nsource
    P(:,source) = repmat(P_mag0(source),K,1).*...
                            exp(1i*(2*pi*(rand(K,1)-0.5)));
end

%___Calculate magnitude and phase of source pressure____
P_mag = abs(P);
P_phi = atan2(imag(P),real(P));

%% Create array data_______________________________________________________
%___Modify PSF mag and/or phase____________________________
PSF_mod = ones(M,nsource);    % Default no modification

%___Construct FFT for each source, data block, and sensor_____________
sourceF = zeros(K,M,nsource);
for source = v_nsource
    for kFFT = v_K
        for m = v_M
            mag = nFFT*P_mag(kFFT,source)*PSF_mod(m,source)/...
                    (8*pi*r_M_nsource(m,source));
            phase = exp(1i*(2*pi*ftrue/c*r_M_nsource(m,source)+...
                            P_phi(kFFT,source)));
            sourceF(kFFT,m,source) = mag.*phase;  
        end
    end
end

%___Sum pressure data for all sources to create array FFT______
arrayF = sum(sourceF,3);

%% Form CSM from FFT data__________________________________________________
%___Calculate strength of incoherent noise___
noise = mean(abs(arrayF(:)))/(10^(SNR/20));

%___Add incoherent noise to array FFT data (uniform distribution used)___
noise_mag = 2*(rand(K,M)-0.5);
noise_phase = 2*(rand(K,M)-0.5);
arrayFn = (noise/mean(abs(noise_mag(:))))*...
            noise_mag.*exp(1i*pi*noise_phase);
arrayF = arrayF + arrayFn;
tK_cTarrayF = two_invK*ctranspose(arrayF);    % Needed for CSC

%___Form G from array FFT data__________
G = two_invK*ctranspose(arrayF)*arrayF;
G_diag = max(0,real(diag(G)));
G(v_G_diag) = G_diag;

%% Create "ambient" CSM needed for CSC processing, which has same 
%   statistics as incoherent noise_________________________________________
amb_mag = 2*(rand(K,M)-0.5);
amb_phase = 2*(rand(K,M)-0.5);
arrayFamb = (noise/mean(abs(amb_mag(:))))*amb_mag.*exp(1i*pi*amb_phase);
Gamb = two_invK*ctranspose(arrayFamb)*arrayFamb;
Gamb_diag = max(0,real(diag(Gamb)));
Gamb(v_G_diag) = Gamb_diag;
aGamb = abs(Gamb); 

%% Form CSM of each source individually without added noise________________
Gs = zeros(M,M,nsource);
for source = v_nsource
    A(:,:) = two_invK*ctranspose(sourceF(:,:,source))*sourceF(:,:,source);
    A(v_G_diag) = max(0,real(diag(A)));
    Gs(:,:,source) = A;     
    clear A;
end


end