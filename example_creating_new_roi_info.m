setenv('FREESURFER_HOME','/Applications/freesurfer');
clc;

if(isempty(which('read_surf','function')))
    addpath([getenv('FREESURFER_HOME') '/matlab']);
end
setenv('FSLDIR','/usr/local/fsl');
if(isempty(which('read_avw','function')))
    addpath([getenv('FSLDIR') '/etc/matlab']);
end

if(isempty(which('spm_vol','function')))
    addpath('/Users/kwj5/MATLAB_TOOLBOXES/spm12');
end
%%


[Avertices, Alabel, Acolortable] = read_annotation([getenv('FREESURFER_HOME') '/subjects/fsaverage/label/lh.aparc.annot']);
[Bvertices, Blabel, Bcolortable] = read_annotation('~/Research/lh.lobes.annot');

lutnames=Acolortable.struct_names;
lutcodes=Acolortable.table(:,5);

fid=fopen('~/Research/FreeSurferROIlist86_Nypipe.txt','r');
M=textscan(fid,'%f=%s');
fclose(fid);
code86=M{1};
name86=M{2};
%%
%this will help with lobe assignment
[A,~,~,~,~,hdr] = read_avw_and_header('/Users/kwj5/Source/nemo/website/atlases/fs86_dil1_allsubj_mode.nii.gz');
A86=A;

%%
%%%%%%

Aidx=find(A(:)>0);
[vi, vj, vk]=ind2sub(size(A),Aidx);
Aijk=[vi vj vk];
Aval=A(Aidx);
uval=unique(Aval);
Axyz=affine_transform(hdr.mat,Aijk);

roixyz=zeros(numel(uval),3);
for i = 1:numel(uval)
    roixyz(i,:)=mean(Axyz(Aval==uval(i),:),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lobe86=cell(size(name86));
hemi86=lobe86;
for i = 1:numel(name86)
    hemi='';
    lobe='';
    lutidx=-1;
    if(regexpimatch(name86{i},'^ctx-lh'))
        hemi='lh';
        [~,lutidx]=ismember(regexprep(name86{i},'^ctx-lh-',''),lutnames);
    elseif(regexpimatch(name86{i},'^ctx-rh'))
        hemi='rh';
        [~,lutidx]=ismember(regexprep(name86{i},'^ctx-rh-',''),lutnames);
    elseif(regexpimatch(name86{i},'^Left'))
        hemi='lh';
        lobe='subcortex';
    elseif(regexpimatch(name86{i},'^Right'))
        hemi='rh';
        lobe='subcortex';
    end
    %{name86{i},lutnames{lutidx} lutcodes(lutidx)}
    %lutcodes(lutidx)
    if(lutidx>0)
        bidx=find(Bcolortable.table(:,5)==mode(Blabel(Alabel==lutcodes(lutidx))));
        lobe=Bcolortable.struct_names{bidx};
    end
    lobe86{i}=sprintf('%s-%s',hemi,lobe);
    hemi86{i}=hemi;
end
name86_nohemi=regexprep(regexprep(name86,'^(ctx-lh-|Left-)',''),'^(ctx-rh-|Right-)','');
name86_nohemi=regexprep(name86_nohemi,'-area|-Proper|-Cortex','');


%make a new version of roixyz points that makes the contralateral pairs exact
%mirror images of each other (cleaner glassbrain)
roixyz_sym=roixyz;
for i = 1:numel(name86)
    if(~strcmp(hemi86{i},'lh'))
        continue;
    end
    otherhemiidx=find(strcmpi(name86_nohemi,name86_nohemi{i}) & strcmpi(hemi86,'rh'));
    if(isempty(otherhemiidx))
        fprintf('No contralateral ROI found for %s\n',name86{i});
        continue;
    elseif(numel(otherhemiidx)>1)
        fprinf('More than one contralateral ROI found for %s. How???\n',name86{i});
        continue;
    end
    newxyz=(roixyz(i,:) + roixyz(otherhemiidx(1),:).*[-1 1 1])/2;
    
    roixyz_sym(i,:) = newxyz;
    roixyz_sym(otherhemiidx(1),:) = newxyz.*[-1 1 1];
end
roixyz=roixyz_sym;

roixyz86=roixyz;

%%

[A,dims,scales,bpp,endian,hdr] = read_avw_and_header('/Users/kwj5/Source/nemo/website/atlases/shen268_MNI1mm_dil1.nii.gz');

%[A,dims,scales,bpp,endian,hdr] = read_avw_and_header('/Users/kwj5/Source/nemo/website/atlases/cc400_new1mm_seq392.nii.gz');

Aidx=find(A(:)>0);
[vi, vj, vk]=ind2sub(size(A),Aidx);
Aijk=[vi vj vk];
Aval=A(Aidx);
uval=unique(Aval);
Axyz=affine_transform(hdr.mat,Aijk);

roixyz=zeros(numel(uval),3);
for i = 1:numel(uval)
    roixyz(i,:)=mean(Axyz(Aval==uval(i),:),1);
end

Aval_to_fs86=zeros(numel(uval),1);
A_to_fs86=zeros(size(A));

lobes=cell(numel(uval),1);
hemis=lobes;
names=lobes;
for i = 1:numel(uval)
    m=mode(A86(A==uval(i) & A86>0));
    isshen268nan=numel(uval)==268 && isnan(m);
    if(isshen268nan)
        %brainstem doesn't math any fs86
        %assign it to cerebellum for now
        return;
        m=18;
        lobes{i}=lobe86{m};
        names{i}=sprintf('%s%d',name86{m},1+sum(regexpmatch(names,['^' name86{m} '[0-9]+$'])));
        fprintf('subcortexnan: %s\n',names{i});
        %lobes{i}='rh-subcortex';
        %names{i}=sprintf('%s%d',name86{m},1+sum(regexpmatch(names,['^' name86{m} '[0-9]+$'])));
        %names{i}='subcortexthing';
    elseif(isnan(m))
        fprintf('isnan; %d\n',i);
    else
        lobes{i}=lobe86{m};
        names{i}=sprintf('%s%d',name86{m},1+sum(regexpmatch(names,['^' name86{m} '[0-9]+$'])));
    end
       
    if(roixyz(i,1)>0)
        hemis{i}='rh';
    else
        hemis{i}='lh';
    end
    Aval_to_fs86(i)=m;
    A_to_fs86(A==uval(i))=Aval_to_fs86(i);
end
names_nohemi=regexprep(regexprep(names,'^(ctx-lh-|Left-)',''),'^(ctx-rh-|Right-)','');
names_nohemi=regexprep(names_nohemi,'-area|-Proper|-Cortex','');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%

lobe_order={'frontal','cingulate','parietal','occipital','temporal','insula','subcortex'};

roiorder={};
roilobe={};
roilobe_idx=[];
roitheta=[];
%thetastep=2*pi/(numel(name86)+2*numel(lobe_order)+4);

thetastep=pi/(sum(regexpmatch(hemis,'lh'))+numel(lobe_order)+2);
theta_prev=pi/2+thetastep;

for i = 1:numel(lobe_order)
    islobe=regexpmatch(lobes,['(lh|Left)-' lobe_order{i}]);
    [~,ysort]=sort(roixyz(islobe,2),'descend');
    yname=names(islobe);
    yname=yname(ysort);
    roiorder=[roiorder; yname];
    roilobe=[roilobe;repmat(lobe_order(i),numel(yname),1)];
    roilobe_idx=[roilobe_idx; repmat(i,numel(yname),1)];
    roitheta=[roitheta; theta_prev+(1:numel(yname))'*thetastep];
    theta_prev=roitheta(end)+thetastep;
end


thetastep=-1*pi/(sum(regexpmatch(hemis,'rh'))+numel(lobe_order)+2);
theta_prev=pi/2+thetastep;
for i = 1:numel(lobe_order)
    islobe=regexpmatch(lobes,['(rh|Right)-' lobe_order{i}]);
    [~,ysort]=sort(roixyz(islobe,2),'descend');
    yname=names(islobe);
    yname=yname(ysort);
    roiorder=[roiorder; yname];
    roilobe=[roilobe;repmat(lobe_order(i),numel(yname),1)];
    roilobe_idx=[roilobe_idx; repmat(i,numel(yname),1)];
    roitheta=[roitheta; theta_prev+(1:numel(yname))'*thetastep];
    theta_prev=roitheta(end)+thetastep;
end


[~,roiorder_idx]=ismember(roiorder,names);
%roilobe=lobe86(roiorder_idx);

roiname_nohemi=regexprep(regexprep(roiorder,'^(ctx-lh-|Left-)',''),'^(ctx-rh-|Right-)','');
roiname_nohemi=regexprep(roiname_nohemi,'-area|-Proper|-Cortex','');
roihemi=regexprep(regexprep(roiorder,'^(ctx-lh|Left).+$','lh'),'^(ctx-rh|Right).+$','rh');

%%
figure;
roixy=[1*cos(roitheta) 1*sin(roitheta)];
plot(roixy(:,1),roixy(:,2),'.k');
islh=regexpmatch(hemis,'rh');
textrad=1.01;
roiname_nohemi=roiorder;
for i = 1:numel(roitheta)
    if(roitheta(i)>=pi/2 && roitheta(i)<=3*pi/2)
        textrotation=roitheta(i)+pi;
        textalign='right';
    else
        textrotation=roitheta(i);
        textalign='left';
    end
    ht=text(textrad*roixy(i,1),textrad*roixy(i,2),roiname_nohemi{i},'rotation',textrotation*180/pi,'horizontalalignment',textalign);
end

[~,roiorder_idx]=ismember(roiorder,names);
%roilobe=lobe86(roiorder_idx);



roixy=[cos(roitheta) sin(roitheta)];
%roixy=[roixy; -roixy(:,1) roixy(:,2)];

roiname_nohemi=regexprep(regexprep(roiorder,'^(ctx-lh-|Left-)',''),'^(ctx-rh-|Right-)','');
roiname_nohemi=regexprep(roiname_nohemi,'-area|-Proper|-Cortex','');
roihemi=regexprep(regexprep(roiorder,'^(ctx-lh|Left).+$','lh'),'^(ctx-rh|Right).+$','rh');

islh=strcmpi(roihemi,'lh');

roilobe_lhrh=strcat(roihemi,'-',roilobe);
roilobe_idx_lhrh=roilobe_idx;
roilobe_idx_lhrh(~islh)=roilobe_idx(~islh)+numel(lobe_order);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%

roi_info=struct('roi_names',{names},'roi_hemi',{hemis},'roi_names_nohemi',{names_nohemi},...
    'roi_xyz',roixyz,'roi_lobes',{lobes},'lobe_order',{lobe_order},...
    'reorder_idx',roiorder_idx,'reorder_roi_theta',roitheta);

save('roi_info_cc400.mat','-struct','roi_info');
