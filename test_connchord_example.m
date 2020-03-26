clear all;
close all;

roi_info=load('roi_info_fs86.mat');
Cstruct=load('hcp_conn_mean.mat');


conn_type='FC';

if(strcmpi(conn_type,'FC'))
    
    C=Cstruct.FCmean;
    
    clim=[-1 1];
    %cmap=GenerateColormap([0 .5 1],[0 0 1; 1 1 1; 1 0 0],256);
    cmap=colormap_roybigbl_gray(256);
elseif(strcmpi(conn_type,'SC'))
    
    C=Cstruct.SCmean;
    
    clim=[0 max(C(:))];
    
    cmap=colormap_roybigbl_gray(256);
    cmap=cmap(129:end,:);
end


fig=connchord(C,roi_info,'clim',clim,'colormap',cmap,'showlabels',true);

%% test glassbrain viewer as well
glassbrain_info=load('glassbrain_surface.mat');

fig=connglass(C,roi_info,glassbrain_info,'clim',clim,'colormap',cmap);

img=connglass(C,roi_info,glassbrain_info,'clim',clim,'colormap',cmap,'threshpercentile',99,'renderviews',{[90 0],[0 90],[0 0]});
figure;
imshow(img);
%imwrite(img,'example_glassbrain_hcp.png');