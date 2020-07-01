function Plot_Results...
            (Y_true,Y_fdbf,Y_csc,Sx,Sy,Sz,pref,sx,sy,RE_switch,...
                SVE_switch,RE_SVE_switch)
            
        
%% Post-process____________________________________________________________
%___Reshape results to match scan grid___
Y_true = reshape(Y_true,Sy,Sx,Sz);
Y_fdbf = reshape(Y_fdbf,Sy,Sx,Sz);
Y_csc = reshape(Y_csc,Sy,Sx,Sz);

%____Convert to dB______________
Y_true = 10*log10(Y_true*pref);
Y_fdbf = 10*log10(Y_fdbf*pref);
Y_csc = 10*log10(Y_csc*pref);

%% Compute variables needed to plot________________________________________
%___Max level used to normalize plots to peak of 0 dB__________
max_lvl = max([max(Y_true(:)),max(Y_fdbf(:)),max(Y_csc(:))]);

%___Colorbar range, dB_____
colorbar_range = [-20 0];  

%___Font size to use for plots____
font_size = 16;

%___Convert inches to cm and specify plot labels_________________________
inch_to_cm = 2.54;
x_tick_labels = round(inch_to_cm*[min(sx(:)),min(sx(:))/2,0,...
                            max(sx(:))/2,max(sx(:))],2);
y_tick_labels = round(inch_to_cm*[max(sy(:)),max(sy(:))/2,0,...
                            min(sy(:))/2,min(sy(:))],2);
x_tick_labels_v = {num2str(x_tick_labels(1)),num2str(x_tick_labels(2)),...
                    num2str(x_tick_labels(3)),num2str(x_tick_labels(4)),...
                    num2str(x_tick_labels(5)),};
y_tick_labels_v = {num2str(y_tick_labels(1)),num2str(y_tick_labels(2)),...
                    num2str(y_tick_labels(3)),num2str(y_tick_labels(4)),...
                    num2str(y_tick_labels(5)),};

%% Plot results____________________________________________________________
%___Plot true source distribution_________________________
subplot(1,3,1);
imagesc(Y_true-max_lvl,colorbar_range);
axis image;
divx = 4;
divy = 4;
box on;  
set(gca,'Xtick',0.5:(Sx-1)/divx:(Sx-0.5),'Xticklabel',...
    x_tick_labels_v,'Ytick',0.5:(Sy-1)/divy:(Sy-0.5),...
    'Yticklabel',y_tick_labels_v,'FontSize',font_size); 
xlabel('X, cm');   
ylabel('Y, cm');
title('Y_T_r_u_e');

%___Plot FDBF results_____________________________________
subplot(1,3,2);
imagesc(Y_fdbf-max_lvl,colorbar_range);
axis image;
box on; 
set(gca,'Xtick',0.5:(Sx-1)/divx:(Sx-0.5),'Xticklabel',...
    x_tick_labels_v,'Ytick',0.5:(Sy-1)/divy:(Sy-0.5),...
    'Yticklabel',{},'FontSize',font_size); 
xlabel('X, cm');   
title('Y_F_D_B_F');

%___Plot CSC2 results______________________________________
subplot(1,3,3);
imagesc(Y_csc-max_lvl,colorbar_range);
axis image;
box on; 
colorbar; 
H = colorbar(); 
set(gca,'Xtick',0.5:(Sx-1)/divx:(Sx-0.5),'Xticklabel',...
    x_tick_labels_v,'Ytick',0.5:(Sy-1)/divy:(Sy-0.5),...
    'Yticklabel',{},'FontSize',font_size); 
xlabel('X, cm');   
if RE_switch == 1
    title('Y_C_S_C (RE constraints)');
elseif SVE_switch == 1
    title('Y_C_S_C (SVE constraints)');
elseif RE_SVE_switch == 1
    title('Y_C_S_C (RE+SVE constraints)');
else
    title('Y_C_S_C (Unconstrained)');
end
set(get(H,'ylabel'),'string','dB','Fontsize',font_size+2);
set(gcf,'Color',[1,1,1]); 


end