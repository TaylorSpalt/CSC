function [S,sx,sy,sz] = ...
            Build_Scan_Grid...
                (sx0,sy0,sz0,Sx,Sy,Sz,gsx,gsy,gsz)


%___Total # of grid points____
S = Sx*Sy*Sz;

%___Transform specified grid into a 3D mesh_______________________
x_start = -gsx*round(Sx/2) + gsx + sx0;
x_end = gsx*round(Sx/2) - gsx + sx0;
y_start = gsy*round(Sy/2) - gsy + sy0;
y_end = -gsy*round(Sy/2) + gsy + sy0;
z_start = -gsz*round(Sz/2) + gsz + sz0;
z_end = gsz*round(Sz/2) - gsz + sz0;
[sx,sy,sz] = ...
    meshgrid...
        (x_start:gsx:x_end,y_start:-gsy:y_end,z_start:gsz:z_end);
 

end