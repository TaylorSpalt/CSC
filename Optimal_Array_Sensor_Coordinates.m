function [mx,my] = ...
            Optimal_Array_Sensor_Coordinates...
                (H,V,rmax,M)

            
% Based on: E.Sarradj, “Optimal planar microphone array arrangements,” 
%               Fortschritte der Akustik-DAGA, 220-223 (2015).

%% Calculate variables needed______________________________________________            
pi_H = pi*H;
v_rmax = 0:0.1:rmax;
l_v_rmax = length(v_rmax);
rmax_base = rmax*sqrt((1:M)/M);
if H >= 0
    p = pi_H*sqrt(1-(v_rmax/rmax).^2);
    pp = pi_H*sqrt(1-(rmax_base/rmax).^2);
else
    p = pi_H*v_rmax/rmax;
    pp = pi_H*rmax_base/rmax;
end

%% Calculate sensor coordinates____________________________________________
z_l_v_rmax = zeros(l_v_rmax,1);
r = zeros(M,1);
for n = 1:M
    sum_m = zeros(n,1);
    for m = 1:n
        Int = z_l_v_rmax;
        for n_int = 1:l_v_rmax
            Int(n_int) = besseli(0,p(n_int));
        end
        sum_m(m) = sum(Int)/besseli(0,pp(m));
    end
    r(n) = sum(sum_m);
end
r = sqrt(r/M);
r = r/max(r)*rmax;
theta = pi*(1:M)'*(1+sqrt(V));
[mx,my] = pol2cart(theta,r);


end