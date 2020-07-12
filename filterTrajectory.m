function [vxyzflt] = filterTrajectory(vxyz, varargin)
% The function filter out the high frequency noise in the trajectory data
% Input variable:
% vxyz: the xyz trajectory, can be n by 3, n by 2, n by 1 or 1 by n array. 

[nflt, tpflt]=init(varargin{:});
isbreak = isnan(vxyz(:,1)); % find where the trajectory are disconnected (by 'nan's)
[labeled, nbreak] = bwlabel(isbreak, 4);

for i=1:nbreak
    indt=find(labeled==i);
    vxyz(indt, :)=vxyz(indt(1)-1, :)+(indt(:)-indt(1)+1)*(vxyz(indt(end)+1, :)-vxyz(indt(1)-1, :))/(indt(end)+1-indt(1)+1);
end
        
vv=diff(vxyz);
ndim=numel(vv)/length(vv);
fflt=zeros(1, nflt);
switch tpflt
    case 'Linear'
        fflt(:)=ones(1, nflt)/nflt;
    case 'Gaussian'
        fflt=fspecial('gaussian', nflt, nflt/5);
        fflt=fflt(1, :)/sum(fflt(1,:));
end

if size(vv, 1)==1; 
    vv=vv';
end

vv=[flipud(vv(2:floor(nflt/2)+1, :)); vv; flipud(vv(end-floor(nflt/2)-1: end-1, :))];
for i=1:ndim
    vv(:,i)=filter(fflt, 1, vv(:,i));
end
vv=vv((floor(nflt/2)+1:end-(floor(nflt/2)))+floor(nflt/2), :);

vxyzflt=cumsum(vv);
vxyzflt=vxyzflt+repmat(mean(vxyz)-mean(vxyzflt), length(vxyzflt), 1);
vxyzflt(isbreak,:)=nan;
end

function [nflt, tpflt]=init(varargin)
for i=2:2:nargin
    switch varargin{i-1}
        case 'NPoint'
            nflt=varargin{i};
        case 'FilterType'
            tpflt=varargin{i};
    end
end
if ~exist('nflt', 'var'), nflt=3; end;
if ~exist('tpflt', 'var'), tpflt='Linear'; end;
end

