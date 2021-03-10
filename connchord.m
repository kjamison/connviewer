function fig=connchord(Cmat,roi_info,varargin)
% fig=connchord(Cmat,roi_info,'param',value,...)
%
% Cmat = NxN matrix or (N*N-N)/2 x 1 upper-triangular vector
%   Note: Cmat must match roi_info in size and ROI order!
% roi_info = struct with ROI display information (see 'roi_info_fs86.mat')
% 
% Params:
% clim: [min max]
% cmap,colormap: eg: jet(256)
% threshval: only display if abs(Cmat)>threshval
% threshpercentile: 0-100. only display if abs(Cmat)>prctile(abs(Cmat), threshpercentile)
% minwidth,maxwidth: line thickness scaled in this range
% showlabels: true (default) = show ROI label around circle
% lobecolor: Lx3 RGB matrix where L=# of lobes (7 for now). default=jet(7)
% axes: if provided, plot on specific axes. otherwise create new figure and return that.
%
% dependencies: regexpmatch, flatten, val2rgb, element

args = inputParser;
args.addParameter('clim',[]);
args.addParameter('cmap',[]);
args.addParameter('colormap',[]);
args.addParameter('lobecolor',[]);
args.addParameter('threshval',[]);
args.addParameter('threshpercentile',90);
args.addParameter('minwidth',0);
args.addParameter('maxwidth',10);
args.addParameter('showlabels',true);
args.addParameter('circlecolor',[.85 .85 .85]);
args.addParameter('backgroundcolor',[1 1 1]);
args.addParameter('axes',[]);

args.parse(varargin{:});
args = args.Results;


numroi=numel(roi_info.roi_names);

name86=roi_info.roi_names;
lobe86=roi_info.roi_lobes;
lobe_order=roi_info.lobe_order;

lobecolor=hsv(numel(lobe_order));
if(isfield(roi_info,'lobe_color') && ~isempty(roi_info.lobe_color))
    lobecolor=roi_info.lobe_color;
end
if(~isempty(args.lobecolor))
    lobecolor=args.lobecolor;
end

minwidth=args.minwidth;
maxwidth=args.maxwidth;

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
if(ischar(cmap))
    cmap=evalin('caller',[cmap '(256)']);
end

showlabels=args.showlabels;
circlecolor=args.circlecolor;
backgroundcolor=args.backgroundcolor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(Cmat,1)==numroi && size(Cmat,2)==numroi)
    %Input should be symmetric
    %If not, check if it is triu or tril
    csym=Cmat-Cmat.';
    
    if(max(abs(csym(:)))>2*eps)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roiorder=roi_info.reorder_idx;
roitheta=reshape(roi_info.reorder_roi_theta,[],1);
roihemi=roi_info.roi_hemi(roiorder);
roiname_nohemi=roi_info.roi_names_nohemi(roiorder);
roixy=[cos(roitheta) sin(roitheta)];
[~,roilobe_idx]=ismember(regexprep(roi_info.roi_lobes(roiorder),'^[rl]h-',''),lobe_order);

islh=strcmpi(roihemi,'lh');


%close all;
if(~isempty(args.axes))
    axes(args.axes);
    fig=args.axes;
else
    fig=figure('color',backgroundcolor);
end

%axcirc=subplot(1,2,2);


textrad=1.025;
hold on;
set(gca,'xlim',[-1 1]*1.25,'ylim',[-1 1]*1.25);

if(showlabels)
    for i = 1:numel(roiorder)
        if(mod(roitheta(i),2*pi)>=pi/2 && mod(roitheta(i),2*pi)<=3*pi/2)
            textrotation=roitheta(i)+pi;
            textalign='right';
        else
            textrotation=roitheta(i);
            textalign='left';
        end
        ht=text(textrad*roixy(i,1),textrad*roixy(i,2),roiname_nohemi{i},'rotation',textrotation*180/pi,'horizontalalignment',textalign);

        set(ht,'color',lobecolor(roilobe_idx(i),:));
        
        %set(ht,'backgroundcolor',[1 1 1]*0);
    end
end
fill(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)),circlecolor,'linestyle','none');
axis off;


% 
% hemis={'lh','rh'};
% for i = 1:numel(lobe_order)
%     for h = 1:numel(hemis)
%         islobe=strcmpi(roilobe_lhrh,[hemis{h} '-' lobe_order{i}]);
%         plot(roixy(islobe,1),roixy(islobe,2),'-','color',lobecolor(i,:),'linewidth',5);
%     end
%     
% end


%text(roixy(islh,1),roixy(islh,2),roiorder(islh),'horizontalalignment','right','color','r');
%text(roixy(~islh,1),roixy(~islh,2),roiorder(~islh),'horizontalalignment','left','color','b');


C_roiorder=Cmat(roiorder,roiorder);

up=triu(ones(size(C_roiorder)),1)>0;

%threshval=0.3;
if(~isempty(args.threshval))
    threshval=args.threshval;
else
    threshval=prctile(abs(C_roiorder(up)),args.threshpercentile);
end
valscale=(abs(C_roiorder)-threshval)/(max(abs(C_roiorder(up)))-threshval);

Cthick=valscale*(maxwidth-minwidth)+minwidth;

Cthresh=abs(C_roiorder)>threshval;

Ccolor=val2rgb(C_roiorder,cmap,clim);

for i = 1:size(C_roiorder,1)
    for j = i+1:size(C_roiorder,2)
        if(~Cthresh(i,j))
            continue;
        end
        p1=roixy(i,:);
        p2=roixy(j,:);

        
        %%%%%%%%%%%%
        if(sqrt(sum((p1-p2).^2))>1.99)
            %if it's basically across the diameter, just draw a line
            path12=[p1; p2];
            path12=[path12; nan nan];
        else
            
            u  = p1';
            v  = p2';
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0
                % ensure the arc is within the unit disk
                theta = [linspace(max(thetaLim),pi,50),...
                    linspace(-pi,min(thetaLim),50)].';
            else
                theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            path12=[r*cos(theta)+x0 r*sin(theta)+y0];
            path12=[path12; nan nan];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %plot(roixy([i j],1),roixy([i j],2),'-','linewidth',Cthick(i,j),'color',flatten(Ccolor(i,j,:)))
        %fill(path12(:,1),path12(:,2),'w','linewidth',Cthick(i,j),'edgecolor',flatten(Ccolor(i,j,:)),'edgealpha',.2);
        fill(path12(:,1),path12(:,2),'w','tag','circpath',...
            'linewidth',Cthick(i,j),'edgecolor',flatten(Ccolor(i,j,:)),'edgealpha',.5*valscale(i,j));
    end
end
plot(roixy(:,1),roixy(:,2),'ko','markerfacecolor','k','tag','circpoint');

set(gca,'clim',clim);
colormap(cmap);
