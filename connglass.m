function figglass = connglass(Cmat,roi_info,glassbrain_info,varargin)
% figglass=connglass(Cmat,roi_info,glassbrain_info,'param',value,...)
%
% Cmat = NxN matrix or (N*N-N)/2 x 1 upper-triangular vector
%   Note: Cmat must match roi_info in size and ROI order!
% roi_info = struct with ROI display information (see 'roi_info_fs86.mat')
% glassbrain_info = struct with glass brain surfaces (see 'glassbrain_surface.mat')
%
% Params:
% view,view_azel: [azimuth elevation] pair in degrees ([90 0]=sagittal, [0 90]=axial, [0 0]=coronal)
% clim: [min max]
% cmap,colormap: eg: jet(256)
% threshval: only display if abs(Cmat)>threshval
% threshpercentile: 0-100. only display if abs(Cmat)>prctile(abs(Cmat), threshpercentile)
% minlinewidth,maxlinewidth: line thickness scaled in this range
% mindotsize,maxdotsize: nodes are scaled in this range (by sum(degree>thresh))
% renderviews: cell array of [azimuth elevation] pairs to render and
% concatenate into a RENDERED return value
%
% Output:
% if 'renderviews' is EMPTY or not specified, return the matlab figure with full 3D interactive GUI
% if 'renderviews' is provided, return rendered images instead
%
% dependencies: flatten, val2rgb, tubeplot, CropBGColor, BuildSphere

args = inputParser;
args.addParameter('view',[]);
args.addParameter('view_azel',[]);
args.addParameter('clim',[]);
args.addParameter('cmap',[]);
args.addParameter('colormap',[]);
args.addParameter('threshval',[]);
args.addParameter('threshpercentile',90);
args.addParameter('minlinewidth',0);
args.addParameter('maxlinewidth',10);
args.addParameter('mindotsize',2);
args.addParameter('maxdotsize',10);
args.addParameter('renderviews',[]);

args.parse(varargin{:});
args = args.Results;

numroi=numel(roi_info.roi_names);

roixyz=roi_info.roi_xyz;


minlinewidth=args.minlinewidth;
maxlinewidth=args.maxlinewidth;

mindotsize=args.mindotsize;
maxdotsize=args.maxdotsize;

clim=[nanmin(Cmat(:)) nanmax(Cmat(:))];
if(~isempty(args.clim))
    clim=args.clim;
end

cmap=jet(256);
if(~isempty(args.cmap))
    cmap=args.cmap;
elseif(~isempty(args.colormap))
    cmap=args.colormap;
end

view_azel=[90 0];
if(~isempty(args.view))
    view_azel=args.view;
elseif(~isempty(args.view_azel))
    view_azel=args.view_azel;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(Cmat,1)==numroi && size(Cmat,2)==numroi)
    %Input should be symmetric
    %If not, check if it is triu or tril
    if(~isequal(Cmat,Cmat.'))
        tmask_up=triu(ones(size(Cmat)),1)>0;
        tmask_down=tmask_up.';
        if(all(Cmat(tmask_down)==0) || all(Cmat(tmask_up)==0))
            Cmat=Cmat+Cmat.';
        else
            error('Input matrix is not symmetric or upper/lower triangular');
        end
    end
elseif(size(Cmat,1)==size(Cmat,2))
    error('Square matrix of size %dx%d does not match roi_info size %d',size(Cmat,1),size(Cmat,2),numroi);
elseif(min(size(Cmat,1),size(Cmat,2))==1)
    Cvec=Cmat;
    n=ceil(sqrt(2*numel(Cvec)));
    if(n ~= numroi)
        error('Input was a vector with length %d, assume square size %dx%d, but roi_info was %d',numel(Cvec),n,n,numroi);
    end

    Cmat=zeros(n,n);
    trimask=triu(ones(n,n),1)>0;
    Cmat(trimask)=Cvec(:);
    Cmat=Cmat+Cmat.';
else
    error('Invalid input size: %dx%d. For roi_info length N, input must be NxN or vector (N*N-N)/2',size(Cmat,1),size(Cmat,2))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isempty(args.renderviews))
    figglass=figure('color',[1 1 1]);
else
    figglass=figure('color',[1 1 1],'visible','off');
end

%subplot(2,2,3);
glassparams={'linestyle','none','facecolor',[1 1 1]*.5,'facealpha',.1};

hpglassL=[];
for i = 1:numel(glassbrain_info.vertices)
    hpglassL=patch('faces',glassbrain_info.faces{i},'vertices',glassbrain_info.vertices{i},glassparams{:});
end

axglass=get(hpglassL(1),'parent');

axis vis3d equal;
axis tight;
hold on;
axis off;

%plot3(roixyz(:,1),roixyz(:,2),roixyz(:,3),'ro');
view(view_azel);


[p,t]=BuildSphere(3);
%diam = 1
fvsphere=struct('faces',t,'vertices',p*.5);
hp=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up=triu(ones(size(Cmat)),1)>0;
if(~isempty(args.threshval))
    threshval=args.threshval;
else
    threshval=prctile(abs(Cmat(up)),args.threshpercentile);
end
valscale=(abs(Cmat)-threshval)/(max(abs(Cmat(up)))-threshval);

Cthick=valscale*(maxlinewidth-minlinewidth)+minlinewidth;

Cthresh=abs(Cmat)>threshval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%roirad=5;

%roirad=nansum(abs(Cmat),2);
roirad=nansum(abs(Cmat.*Cthresh),2);
roirad_mask=~isnan(roirad) & roirad>0;

roirad=(roirad-min(roirad))/(max(roirad)-min(roirad));

roirad=roirad*(maxdotsize-mindotsize)+mindotsize;



for i = 1:size(roixyz,1)
    if(~roirad_mask(i))
        continue;
    end
    sphvert=fvsphere.vertices;
    if(numel(roirad)>1)
        sphvert=sphvert*roirad(i);
    else
        sphvert=sphvert*roirad;
    end
    
    sphvert=bsxfun(@plus,sphvert,roixyz(i,:));
    hp(i)=patch('faces',fvsphere.faces,'vertices',sphvert,...
        'linestyle','none','facecolor','k','facealpha',.5);
end
lighting phong;
material dull;


Ccolor=val2rgb(Cmat,cmap,clim);

for i = 1:size(Cmat,1)
    for j = i+1:size(Cmat,2)
        if(~Cthresh(i,j))
            continue;
        end
        p1=roixyz(i,:);
        p2=roixyz(j,:);
        path12=[p1; mean([p1; p2],1); p2];
        %fill3(path12(:,1),path12(:,2),path12(:,3),'w','linewidth',Cthick(i,j)+2,'edgecolor',flatten(Ccolor(i,j,:)),'edgealpha',.5*valscale(i,j));
        [hcyl, hsph] = tubeplot(path12,(Cthick(i,j)+2)/4,[0 0 0],{'spheres',false},...
            'facecolor',flatten(Ccolor(i,j,:)),'linestyle','none','facealpha',.5*valscale(i,j));
    end
end
view(view_azel);

set(gca,'clim',clim);
colormap(cmap);


if(~isempty(args.renderviews))
    allviews=args.renderviews;
    
    tmpd=tempname;
    mkdir(tmpd);
    
    
    
    fig_offscreen=figure('visible','off','position',get(0,'screensize'),'color',[1 1 1]);
    newax=copyobj(axglass,fig_offscreen);
    
    %close(figglass);
    
    img_allview={};
    for v = 1:numel(allviews)
        
        view(newax,allviews{v});
        drawnow;
        
        tmpimgfile=sprintf('%s/glassbrain_%d.png',tmpd,v);
        saveas(fig_offscreen,tmpimgfile);
        img_allview{end+1}=imread(tmpimgfile);
    end
    
    rmdir(tmpd,'s');
    close(fig_offscreen);
    
    croprect={};
    for i = 1:numel(img_allview)
        [~,croprect{i}] = CropBGColor(img_allview{i},img_allview{1}(1,1,:));
    end
    croprect=cat(1,croprect{:});
    croprect=[min(croprect(:,1:2)) max(croprect(:,3:4))];
    for i = 1:numel(img_allview)
        img_allview{i}=img_allview{i}(croprect(1):croprect(3),croprect(2):croprect(4),:);
    end
    
    img_new=cat(2,img_allview{:});

    close(figglass);

    figglass=img_new;
end

