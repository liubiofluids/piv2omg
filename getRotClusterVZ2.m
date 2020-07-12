function [Omega, z, uout, vout, anglexy]=getRotClusterVZ2(x, y, u, v, varargin)
[cntrrot, errtol, vstep, stepdiff, zdep, omega_z, ang0, mskrt] = init(x, y, varargin{:});
n=numel(x(:));
if isempty(cntrrot)
    cntrrot=[mean(x(~isnan(u))), mean(y(~isnan(u))), sqrt(std(x(~isnan(u)))^2+std(y(~isnan(u)))^2)];
end
x0=cntrrot(1);
y0=cntrrot(2);
dz=cntrrot(3);
[Omega, z, rang] = getOmega(x, y, u, v, x0, y0, zdep, 'Omega_z', omega_z, 'AngleXY', ang0, 'MaskRatio', mskrt);

[uout, vout] = getuv(x, y, z, Omega, x0, y0, 0);
uout=reshape(uout, size(x));
vout=reshape(vout, size(x));
anglexy=rang;
end



function [cntrrot, errtol, vstep, stepdiff, zdep,  omega_z, ang0, mskrt] = init(x, y, varargin)
n=nargin;
for i=2:2:n-1
switch varargin{i-1}
case 'GuessCenter'
cntrrot=varargin{i}; 
case 'Tolerance'
errtol=varargin{i};
case 'Step'
vstep=varargin{i};
case 'StepDifference'
stepdiff=varargin{i};
case 'ZDepth'
zdep=varargin{i};
case 'Omega_z'
omega_z=varargin{i};
case 'AngleXY'
        ang0=varargin{i};
case 'MaskRatio'
mskrt=varargin{i}; 

end
end


if ~exist('cntrrot', 'var'), cntrrot=[]; end;
if ~exist('errtol', 'var'), errtol=1E-10; end;
if ~exist('vstep', 'var'), vstep=ones(1,3)*1E-2; end;
if ~exist('stepdiff', 'var'), stepdiff=1; end;
if ~exist('zdep', 'var'), zdep=10; end;
if ~exist('omega_z', 'var'), omega_z=[]; end;
if ~exist('ang0', 'var'), ang0=0; end;
if ~exist('mskrt', 'var'), mskrt=1; end;
end

function [Omega, z, rang] = getOmega(x, y, u, v, x0, y0, dz, varargin)

[omega_z, ang0, mskrt] = initsub(varargin{:});
vr=sqrt((x(:)-mean(x(:))).^2+(y(:)-mean(y(:))).^2);
[vrsrt, si]=sort(vr);
nanmsk=ones(size(x));
nanmsk(si(floor(mskrt*length(si))+1:length(si)))=nan;
u1=u; v1=v; x1=x; y1=y;
x(isnan(u))=nan;
y(isnan(u))=nan;
v=v.*nanmsk;
u=u.*nanmsk;
x=x.*nanmsk;
y=y.*nanmsk;
vm=nanmean(v(:));
um=nanmean(u(:));
xm=nanmean(x(:));
ym=nanmean(y(:));
Omega=zeros(3,1);
if isempty(omega_z)
Omega(3)=(nanmean(v(:).*x(:))-nanmean(u(:).*y(:))-vm*xm+um*ym)/(nanmean(x(:).^2)-xm^2+nanmean(y(:).^2)-ym^2);
else
Omega(3)=omega_z;
end
uex=u+(y-y0)*Omega(3);
vex=v-(x-x0)*Omega(3);

%quiver(x(:), y(:), uex(:), vex(:));return;
func=@(x) norm(evalEqu(x, uex, vex)); %scoring function
[rang] = fminsearch(func, ang0);
%z=.5*(uex/sin(rang)-vex/cos(rang));
uex1=u1+(y1-y0)*Omega(3);
vex1=v1-(x1-x0)*Omega(3);
z=(uex1*(sin(rang))-vex1*cos(rang));
Omega(1:2)=[cos(rang), sin(rang)]*nanstd(z(:))/dz;
z=z*dz/nanstd(z(:));

function [omega_z, ang0, mskrt]=initsub(varargin)

for i=2:2:nargin
switch varargin{i-1}
case 'Omega_z'
omega_z=varargin{i};
    case 'AngleXY'
        angle0=varargin{i};
case 'MaskRatio'
	mskrt=varargin{i}; %% including only center u information for Omega and z 
end
end

if ~exist('omega_z', 'var'), omega_z=[]; end;
if ~exist('ang0', 'var'), ang0=0; end;
if ~exist('mskrt', 'var'), mskrt=1; end;

end
end


function [u, v] = getuv(x, y, z, Omega, x0, y0, dz)
n=numel(x(:));
vv3d=cross(repmat(Omega(:)', n, 1),[x(:)-x0, y(:)-y0, z(:)-dz*ones(n, 1)], 2);
u=vv3d(:,1);
v=vv3d(:,2);
end

function [verr] = evalEqu(angxy, uex, vex)
msin=repmat(sin(angxy(:))',  numel(uex), 1);
mcos=repmat(cos(angxy(:))', numel(vex), 1);
mcos_sin2=repmat((cos(angxy(:))./sin(angxy(:)).^2)',  numel(uex), 1);
msin_cos2=repmat((sin(angxy(:))./cos(angxy(:)).^2)',  numel(vex), 1);
mu=repmat(uex(:), 1, numel(angxy));
mv=repmat(vex(:), 1, numel(angxy));
verr=[nansum((mu.*mcos+mv.*msin).^2)]; %./(mu.^2+mv.^2))];
end


