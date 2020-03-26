function [hcylinder hsphere] = tubeplot(verts,widths,cvalues,args,varargin)


%rotate_sphere = false;

%%%%%%%%%%%%
parser = inputParser;
parser.addParameter('resolution',20);
parser.addParameter('interp',1);
parser.addParameter('spheres',true);
parser.addParameter('cylinders',true);
parser.addParameter('spherescale',.97);

parser.parse(args{:});
argstruct = parser.Results;

cylres=argstruct.resolution;
interpfactor=argstruct.interp;
show_sphere=argstruct.spheres;
show_cyl=argstruct.cylinders;
spherescale=argstruct.spherescale;

%%%%%%%%%%%%


if(nargin < 3 || isempty(cvalues))
    cvalues = [];
end
cvalues = cvalues(:); %make sure it's a column vector


N = size(verts,1);

if(numel(widths) == 1)
    widths = widths*ones(N,1);
end

verts = spline(linspace(0,1,N),verts',linspace(0,1,interpfactor*N))';
widths = spline(linspace(0,1,N),widths',linspace(0,1,interpfactor*N));
cvalues = spline(linspace(0,1,N),cvalues',linspace(0,1,interpfactor*N));

D = diff(verts);
Dlen = sqrt(sum(diff(verts).^2,2));

horz_idx = D(:,3) == 0;

Daxis = cross(D,[D(:,1:2) zeros(size(D,1),1)],2);
Daxis_horz = cross(D,[D(:,1:2) -ones(size(D,1),1)],2);

Daxis(horz_idx,:) = Daxis_horz(horz_idx,:); %need other crossprod reference for horz lines
Daxis_len = sqrt(sum(Daxis.^2,2));

%Daxis(Daxis_len == 0,:) = repmat([1 0 0],sum(Daxis_len == 0), 1);
Daxis_len(Daxis_len == 0) = 1;
Daxis = Daxis./repmat(Daxis_len,[1 3]);


Dtheta = sign(D(:,3)).*acos(D(:,3)./Dlen);
Dtheta(horz_idx) = acos(0); %horz lines are always pi;
%%
hsph = [];
hcyl = [];
if(show_sphere)
    
    hsph = plotsphere(verts,widths*spherescale,cvalues,'resolution',cylres,varargin{:});
end

%[Xs Ys Zs] = sphere(cylres);
[X Y Z] = cylinder(1,cylres);
%X = X-mean(X(:));
hold on;

if(show_sphere)
    hsph = zeros(size(verts,1),1);
end

if(show_cyl)
    hcyl = zeros(size(verts,1)-1,1);
end

for i = 1:size(verts,1)-1
    xi = diag(widths(i:i+1))*X;
    yi = diag(widths(i:i+1))*Y;
    zi = diag([0 Dlen(i)])*Z;
    
    %xsi = spherescale*width(i)*Xs;
    %ysi = spherescale*width(i)*Ys;
    %zsi = spherescale*width(i)*Zs;

    [ptmp R] = RodrigRotate(Daxis(i,:),Dtheta(i),[0 0 1]);
    
    p = [R*[xi(:) yi(:) zi(:)]']';
    p = p + repmat(verts(i,:),size(p,1),1);
    xi = reshape(p(:,1),2,size(p,1)/2);
    yi = reshape(p(:,2),2,size(p,1)/2);
    zi = reshape(p(:,3),2,size(p,1)/2);
    
    %ps = [R*[xsi(:) ysi(:) zsi(:)]']';
    %ps = ps + repmat(verts(i,:),size(ps,1),1);
    %xsi = reshape(ps(:,1),size(xsi));
    %ysi = reshape(ps(:,2),size(ysi));
    %zsi = reshape(ps(:,3),size(zsi));
    
    if(show_cyl)
        hcyl(i) = surface(xi,yi,zi,varargin{:});
        if(~isempty(cvalues))
            set(hcyl(i),'CData',repmat(cvalues(i:i+1)',1,cylres+1));
        end
    end
    %hsph(i) = surface(xsi,ysi,zsi,'facecolor','r','facealpha',1,'edgealpha',0);
end

if(nargout > 0)
    hcylinder = hcyl;
end

if(nargout > 1)
    hsphere = hsph;
end


function [pnew Rmat] = RodrigRotate(a,theta,p)

A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
R = eye(3) + A*sin(theta) + A*A*(1-cos(theta));
R = expm(A*theta);
%k = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]
%R = eye(3) + k*sin(theta) + (1-cos(theta))*k*k

if(nargin > 2 && ~isempty(p))
    pnew = [R*p']';
else
    pnew = [];
end

if(nargout > 1)
    Rmat = R;
end

%%
function hsurf = plotsphere(pos,radius,cvalues,varargin)

if(nargin < 3 || isempty(cvalues))
    cvalues = [];
end

if(numel(varargin) > 1 && strcmpi(varargin{1},'resolution'))
    sphereres = varargin{2};
    varargin = {varargin{3:end}};
else
    sphereres = 20;
end

if(min(size(pos)) == 1)
	pos = reshape(pos,1,3);
end

[sx sy sz] = sphere(sphereres);
np = get(gca,'nextplot');
set(gca,'nextplot','add');
hsurf = zeros(size(pos,1),1);
for i = 1:size(pos,1)
    nr = min(i,numel(radius));
    hsurf(i) = surface(sx*radius(nr)+pos(i,1), sy*radius(nr)+pos(i,2), sz*radius(nr)+pos(i,3),varargin{:});
    if(~isempty(cvalues))
		nc = min(i,size(cvalues,1));
		if(size(cvalues,2) == 1)
			cv = repmat(cvalues(nc),size(get(hsurf(i),'CData')));
		elseif(size(cvalues,2) == 3)
			cv = repmat(reshape(cvalues(nc,:),[1 1 3]),size(get(hsurf(i),'CData')));
		end
        set(hsurf(i),'CData',cv);
    end
end
set(gca,'nextplot',np);