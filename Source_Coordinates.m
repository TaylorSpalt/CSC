function [xsource,ysource,zsource,nsource] = ...
            Source_Coordinates...
                (spert,R)


%% Define source coordinates in Cartesian system___________________________
xsource = [0;7;1;-4];        	% inches      
ysource = [-1;2;-7;5];      	% inches
nsource = length(xsource);  
zsource = R*ones(nsource,1);	% inches

%% Randomly perturb source positions if selected (Uniform distribution)____
if spert > 0
    xsource = xsource + (rand(nsource,1)-0.5)*spert;
    ysource = ysource + (rand(nsource,1)-0.5)*spert;
    zsource = zsource + (rand(nsource,1)-0.5)*spert;
end


end