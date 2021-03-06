function [varargout] = Cluster2PIV(folder, flname, flcord, fldisp, winsz, varargin)
%CLUSTER2PIV PIV analyses and associated rigid-body rotation
%   [vomg, vz] = Cluster2PIV(folder, flname, flcord, fldisp, winsz, OPTIONS);
%   computes the angular velocity and axial variance of imaged particles 
%   associated with a quasi-2D velocity field due to the finite focal depth
%   in imaging. 
%
%   Prerequisite Packages
%   ---------------------
%   PIVlab: Available from MathWorks File Exchange
%
%   Class Support
%   -------------
%   folder, flname, flcord, and fldisp are strings.
%       flcord and fldisp are the outputs from running the 
%       'trackbackground' GUI. 
%   winsz = [m, n] is a 1x2 array of integers.
%   vomg is a k x 3 double array.
%   vz = is a 3D double array, e.g., k x m x n. 
%   r2zdepthdat.mat in the working folder is the specific calibration data,
%   which should be customized to the corresponding imaging setup.
%
%   Example
%   -------
%   [vomg] = Cluster2PIV('./testframes/', 'testframe%d.jpg',
%   './testdata/vcord.txt', './testdata/vdisp.txt', [120, 120]);
%   
%   

varargout{1}=[];
[fnmov, xyscale, btime, t0, vframe, vcrop, vresize, nflt, tpflt, rmask, pivconf, outfl, iquad, pivfl, bomg, zdep, npres, mskrt, mskfl, brel, bomgz, outflfit, stpvaromg, indsel] = init(varargin{:});

if ~isempty(fnmov)
aviobj=VideoWriter(fnmov); aviobj.FrameRate=12;
aviobj.open();
end

if ~isnumeric(flcord) && ~isnumeric(fldisp)
vcord=cleardup(load(flcord));
vdisp=cleardup(load(fldisp));
vdisp0=vdisp;
vcord0=vcord;

if nflt>0
vpos=filterTrajectory(vcord(:, 2:3)-cumsum(vdisp(:,2:3)), 'NPoint', nflt, 'TypeFilter', tpflt);
vcord(:,2:3)=cumsum(vdisp(:,2:3))+vpos;
vcordsm=vcord; 
vcordsm(:, 2:3)=filterTrajectory(vcord(:, 2:3), 'NPoint', nflt, 'TypeFilter', tpflt);
else
vcordsm=vcord;
end
else
vind=flcord;
flcord=[fldisp(1:max(strfind(fldisp, '/'))), 'vcord.txt'];
vcord=cleardup(load(flcord));
vdisp=cleardup(load(fldisp));
vpos=filterTrajectory(vcord(:, 2:3)-cumsum(vdisp(:,2:3)), 'NPoint', nflt, 'TypeFilter', tpflt);
vcord(:,2:3)=cumsum(vdisp(:,2:3))+vpos;
vcordsm=vcord; 
vcordsm(:, 2:3)=filterTrajectory(vcord(:, 2:3), 'NPoint', nflt, 'TypeFilter', tpflt);
vcord=vcord(vcord(:,1)>=min(vind) & vcord(:,1)<=max(vind), :);
vcordsm=vcordsm(vcordsm(:,1)>=min(vind) & vcordsm(:,1)<=max(vind), :);
end

if rmask>0
imgmask=zeros(fliplr(winsz));
[xm, ym]=meshgrid(1:winsz(1), 1:winsz(2));
xm=xm-mean(xm(:));
ym=ym-mean(ym(:));
rm=sqrt(xm.^2+ym.^2);
imgmask(rm<rmask)=1;
imgmask(floor(size(rm, 1)/2+1), :)=1; % to avoid a simply connected mask
[B,L] = bwboundaries(1-imgmask,'holes');
for k=1:length(B)
    pivconf.mask_inpt{k,1}=B{k}(:,1);
    pivconf.mask_inpt{k,2}=B{k}(:,2);
end
n1=length(B);
imgmask(:)=0;
imgmask(rm<rmask)=1;
imgmask(floor(size(rm, 1)/2-1), :)=1;
[B,L] = bwboundaries(1-imgmask,'holes');
for k=1:length(B)
    pivconf.mask_inpt{k+n1,1}=B{k}(:,1);
    pivconf.mask_inpt{k+n1,2}=B{k}(:,2);
end
end

img=imread([folder, sprintf(flname, vcord(1, 1))]);

set(gcf, 'Unit', 'pixel', 'Position', [100, 100, size(img, 2), size(img, 1)]);
set(gca, 'Position', [0, 0, 1, 1]);
if fliplr(size(img))==winsz
    rectROI=[1,1, winsz];
else
    rectROI=[vcord(1, 2:3)-.5*winsz, winsz];
end
imgcrp1=PIVlab_preproc(imcrop(img, rectROI), [], 0, 30, 0, 15,...
    0, 0, 3, [], []);
if ~isempty(mskfl)
    imgmask=imread(sprintf(mskfl, vcord(1, 1)));
    imgmaskcrp1=imcrop(imgmask, rectROI);
end

if btime
    data=load([folder, '../data_cord.txt']);
    if t0==0
        if size(vcord, 1)>2
            if abs(diff(data(vcord(1:2, 1), 1)))>10*abs(diff(data(vcord(2:3, 1), 1)))
                t0=data(vcord(2, 1),1);
            
            else
                t0=data(vcord(1,1), 1);
            end
        else
            t0=data(vcord(1,1), 1);
        end
    end
    vt=t0-data(vcord(:,1), 1);
end

if isempty(indsel)
    ilist=1:length(vcord(:,1))-1;
else
    imin=find(vcord(:,1)==min(indsel));
    imax=find(vcord(:,1)==max(indsel));
    if isempty(imin), imin=1; end;
    if isempty(imax), imax=length(vcord(:,1))-1; end;
    ilist=imin:imax;
end

for i=ilist
img=imread([folder, sprintf(flname, vcord(i+1, 1))]);
if fliplr(size(img))==winsz
    rectROI=[1,1, winsz];
else
    rectROI=[vcord(i, 2:3)-.5*winsz, winsz];
end
imgcrp2=PIVlab_preproc(imcrop(img, rectROI), [], 0, 30, 0, 15,...
    0, 0, 3, [], []);
if norm(size(imgcrp2)-size(imgcrp1))>0
    img=imread([folder, sprintf(flname, vcord(i, 1))]);
    imgcrp1=PIVlab_preproc(imcrop(img, rectROI), [], 0, 30, 0, 15,...
    0, 0, 3, [], []);    
end

if isempty(pivfl)
    if ~isempty(mskfl)
        imgmask=imread(sprintf(mskfl, vcord(i, 1)));
        imgmaskcrp2=imcrop(imgmask, rectROI);
        if norm(size(imgmaskcrp2)-size(imgmaskcrp1))>0
            imgmask=imread(sprintf(mskfl, vcord(i-1, 1)));
            imgmaskcrp1=imcrop(imgmask, rectROI);
        end
        if brel
            [xg, yg]=meshgrid(1:size(imgmaskcrp1, 2), 1:size(imgmaskcrp1, 1));
            COM1=[mean(xg(imgmaskcrp1(:)>0)), mean(yg(imgmaskcrp1(:)>0))];
            COM2=[mean(xg(imgmaskcrp2(:)>0)), mean(yg(imgmaskcrp2(:)>0))];
            vshft=COM1-COM2;
            imgmaskcrp1=imshftp(double(imgmaskcrp1), vshft)>0;
            imgcrp1=imshftp(double(imgcrp1), vshft);
        end
        imgmaskf=imgmaskcrp2&imgmaskcrp1;
        
        [labeled, nmsk]=bwlabel(imgmaskf, 4);
        graindata=regionprops(labeled);
        vym=zeros(1, nmsk);
        for ii=1:nmsk
            vym(ii)=graindata(ii).Centroid(2);
        end
        pivconf.mask_inpt={};
        imgmaskf(max([floor(vym-1); ones(1, nmsk)]), :)=1;
        [B,L] = bwboundaries(1-imgmaskf,'holes');
        for k=1:length(B)
            pivconf.mask_inpt{k,1}=B{k}(:,2);
            pivconf.mask_inpt{k,2}=B{k}(:,1);
        end
        imgmaskf=imgmaskcrp2&imgmaskcrp1;
        imgmaskf(min([floor(vym+1); ones(1, nmsk)*winsz(2)]), :)=1;        
        n1=length(B);
        [B,L] = bwboundaries(1-imgmaskf,'holes');
        for k=1:length(B)
            pivconf.mask_inpt{n1+k,1}=B{k}(:,2);
            pivconf.mask_inpt{n1+k,2}=B{k}(:,1);
        end
    end
   
[x, y, u, v, typevector] = piv_FFTmulti (imgcrp1,imgcrp2,pivconf.interrogationarea, pivconf.step, ...
    pivconf.subpixfinder, pivconf.mask_inpt, pivconf.roi_inpt, pivconf.passes,...
    pivconf.int2, pivconf.int3, pivconf.int4, pivconf.imdeform, pivconf.repeat, ...
    pivconf.mask_auto, 0);
fprintf('\b\b'); % stop displaying too much informaion
else
    load(sprintf(pivfl, vcord(i, 1)));
end
imgcrp1=imgcrp2;
if exist('imgmaskcrp2', 'var')
imgmaskcrp1=imgmaskcrp2;
end

pfull=axes('parent', gcf, 'position', [0, 0, 1, 1]);
himg1=imshow(adapthisteq(mat2gray(img)));
if rmask>0
hcirc=viscircles(vcord(i, 2:3), rmask, 'LineWidth', 1, 'LineStyle', '--');
end

if isempty(iquad)
    iinset= (vcordsm(i,2)<size(img, 2)/2)*2+vcordsm(i, 3)<size(img, 1)/2;
else
    iinset=iquad;
end

switch(iinset)
    case 4
        pinst=axes('parent', gcf, 'position', [0, 0, 1, 1]);
    case 3
        pinst=axes('parent', gcf, 'position', [0.6, 0, 0.4, 0.4*size(img, 2)/size(img, 1)]);
    case 2
        pinst=axes('parent', gcf, 'position', [0.6, 1-0.4*size(img, 2)/size(img, 1), 0.4, 0.4*size(img, 2)/size(img, 1)]);
    case 1
        pinst=axes('parent', gcf, 'position', [0, 0, 0.4, 0.4*size(img, 2)/size(img, 1)]);
    case 0
        pinst=axes('parent', gcf, 'position', [0, 1-0.4*size(img, 2)/size(img, 1), 0.4, 0.4*size(img, 2)/size(img, 1)]);
end

himg2=imshow(adapthisteq(mat2gray(imgcrp2)));
if bomg
    if ~exist('mz', 'var')
        mz=zeros(length(vcord(:,1))-1, size(u, 1), size(u, 2));
    end
    if isempty(stpvaromg)
        [Omega, z, uout, vout, anglexy] = getRotClusterVZ2(x, y, u, v, 'ZDepth', r2zdepth(rmask*xyscale)/xyscale, 'NPresearch', npres, 'MaskRatio', mskrt);
        vOmega(i, :)=Omega;
        mz(i, 1:size(z,1), 1:size(z,2))=z;
    else
        [ix0, iy0] = meshgrid(1:size(x, 2), 1:size(x, 1));
        [ix, iy] = meshgrid(1:stpvaromg(2):size(x, 2), 1:stpvaromg(2):size(x, 1));
        spd=sqrt(u.^2+v.^2);
        spd(isnan(spd))=0;
        [labeled, nspd] = bwlabel(spd>0, 4);
        vpospole=[];
 
        for ii=1:nspd
            spdblr=im2blur(spd.*(labeled==ii), 5); spdblr=(spdblr-min(spdblr(:)))/(max(spdblr(:))-min(spdblr(:)));
            bwspd=im2bw(spdblr, 0.9);
            [labeledp, npl] = bwlabel(bwspd, 4);
            for jj=1:npl
                vpospole=[vpospole; [mean(x(labeledp==jj)), mean(y(labeledp==jj))]];
            end
        end
        mures=[];
        mvres=[];
        muout=[];
        mvout=[];
        mvaromg=[];
        mvarz=[];
        vanglexy=[];
        for ii = 1:numel(ix(:))
            ROI=(ix0-ix(ii))>=0 & (ix0-ix(ii))<stpvaromg(1) & (iy0-iy(ii))>=0 & (iy0-iy(ii))<stpvaromg(1);
            uroi=u; vroi=v;
            uroi(~ROI)=nan;
            vroi(~ROI)=nan;
            if sum(~isnan(u(ROI)))>stpvaromg(3)
                xm=mean(x(~isnan(uroi)));
                ym=mean(y(~isnan(uroi)));
                vdist=sum((vpospole-repmat([xm, ym], size(vpospole, 1), 1)).^2, 2);
                [sd, si]=sort(vdist);
                gusscntr=[vpospole(si(1), :), sqrt(std(x(~isnan(uroi)))^2+std(y(~isnan(uroi)))^2)];
                [Omega, z, uout, vout, anglexy] = getRotClusterVZ2(x, y, uroi, vroi, 'ZDepth', r2zdepth(rmask*xyscale)/xyscale, 'NPresearch', npres, 'GuessCenter', gusscntr);
                mvaromg=[mvaromg; Omega];
                mvarz=cat(3, mvarz, z);
                muout=cat(3, muout, uout);
                mvout=cat(3, mvout, vout);
                vanglexy=[vanglexy; anglexy];
                mures=cat(3, mures, Omega(2)*z);
                mvres=cat(3, mvres, -Omega(1)*z);
            end
        end
        z=nanmean(mvarz, 3);
        uout=nanmean(muout, 3);
        vout=nanmean(mvout, 3);
        anglexy=nanmean(anglexy);
        ures=nanmean(mures, 3);
        vres=nanmean(mvres, 3);
        Omega=nanmean(mvaromg);
        stdOmega=sqrt(sum(nanmean((mvaromg-repmat(Omega, size(mvaromg, 1), 1)).^2)));
        vOmega(i, :)=Omega; 
    end
end

if bomgz && exist('Omega', 'var')
    if isempty(stpvaromg)
        u=u-Omega(2)*z;
        v=v+Omega(1)*z;
        uout=uout-Omega(2)*z;
        vout=vout+Omega(1)*z;
    else
        u=u-ures;
        v=v-vres;
        uout=uout-ures;
        vout=vout-vres;
    end
end

hold on; hq=quiver(x, y, 10*u, 10*v); 
set(hq, 'AutoScale', 'on', 'Color', [1, .2, .2], 'LineWidth', 2, 'AutoScaleFactor', 1.6);
if bomg
	hqf=quiver(x, y, 10*uout, 10*vout); 
    set(hqf, 'AutoScale', 'on', 'Color', [.2, .2, 1], 'LineWidth', 2, 'AutoScaleFactor', 1.6);
end
hrct=rectangle('Position',[1, 1, size(imgcrp2, 2)-1, size(imgcrp2, 1)-1],'FaceColor', 'none','EdgeColor', 'w');
if btime
tb=annotation('textbox', 'Units', 'Normalize', ...
    'Position',[0.06, 0.92, 0.2, 0.06*size(img, 1)/size(img, 2)], 'BackgroundColor', 'k', ...
    'FontSize', 18, 'FontWeight', 'Bold', 'Color', [1, 1, 1], ...
    'VerticalAlignment', 'middle', 'Margin', 10);
tb.String=sprintf('t = %.2f s', vt(i));
%t=text(0.1, 0.95, sprintf('t = %.2f s', vt(i)), 'Units', 'Normalized', 'FontSize', 20, 'FontWeight', 'Bold', 'Color', [1, 1, 1]);
set(gca, 'Color', 'k');
end
figure(gcf); pause(0.1);
if exist('aviobj', 'var')
hmv=getframe(gcf);
frm=hmv.cdata;
if ~isempty(vcrop)
    frm=imcrop(frm, vcrop);
end
if ~isempty(vresize)
    frm=imresize(frm, vresize);
end
writeVideo(aviobj, frm);
end
if i<length(vcord(:,1))-1
delete(hq);
delete(himg1);
delete(himg2);
delete(hrct);
if exist('hcirc', 'var'), delete(hcirc); end; 
if exist('tb', 'var'), delete(tb); end;
if exist('hqf', 'var'), delete(hqf); end;
delete(pfull);
delete(pinst);
end
if ~isempty(outfl)
    save(sprintf(outfl, vcord(i, 1)), 'x', 'y', 'u', 'v', 'typevector');
    if ~isempty(stpvaromg)
        save(sprintf(outfl, vcord(i, 1)), 'stdOmega', '-append');
    end
end
if ~isempty(outflfit) && exist('uout', 'var')
    save(sprintf(outflfit, vcord(i, 1)), 'x', 'y', 'uout', 'vout', 'typevector');
end
hold off;
end

%imagesc(imgcrp2); hold on;
%quiver(x, y, u, v); return;
if bomg
    varargout{1}=vOmega;
    varargout{2}=mz;
end

if exist('aviobj', 'var')
close(aviobj);
end

end

function [fnmov, xyscale, btime, t0, vframe, vcrop, vresize, nflt, tpflt, rmask, pivconf, outfl, iquad, pivfl, bomg, zdep, npres, mskrt, mskfl, brel, bomgz, outflfit, stpvaromg, indsel] = init(varargin)
pivconf=struct(...
    'interrogationarea', 16, ...
    'step', 16, ...
    'subpixfinder', 1, ...
    'mask_inpt', [], ...
    'roi_inpt', [],...
    'passes', 2, ...
    'int2', 16, ...
    'int3', 16, ...
    'int4', 16, ...
    'imdeform', '*spline',...
    'repeat', 0,...
    'mask_auto', 0);
for i=2:2:nargin
switch varargin{i-1}
case 'MovieFile'
fnmov = varargin{i};
    case 'XYScale'
        xyscale=varargin{i};
    case 'ShowTime'
        btime=varargin{i};
    case 'Time0'
        t0=varargin{i};
    case 'SelectFrame'
        vframe=varargin{i};
    case 'Resize'
        vresize=varargin{i};
    case 'Crop'
        vcrop=varargin{i};
    case 'NPoint'
        nflt=varargin{i};
    case 'FilterType'
        tpflt=varargin{i};
    case 'MaskRadius'
        rmask=varargin{i};
    case 'OutFile'
        outfl=varargin{i};
    case 'InsetQuad'
        iquad=varargin{i};
    case 'PIVData'
        pivfl=varargin{i};
    case 'FitOmega'
        bomg=varargin{i};
    case 'ZDepth'
        zdep=varargin{i};
    case 'MaskRatio' % include center data for PIV
        mskrt=varargin{i};
    case 'MaskFile'
        mskfl=varargin{i};
    case 'Relative'
        brel=varargin{i}; % calculate the relative velocity in addition to COM
    case 'Zcomponent'
        bomgz = varargin{i}; % find the velocity field due to omgz
    case 'OutFileFit'
        outflfit=varargin{i};
    case 'VarOmega'
        stpvaromg=varargin{i};
    case 'SelectIndex'
        indsel=varargin{i};
end
end
if ~exist('fnmov', 'var'), fnmov=[]; end;
if ~exist('radblur', 'var'), radblur=20; end;
if ~exist('wghtmax', 'var'), wghtmax=1; end;
if ~exist('thrshld', 'var'), thrshld=.1; end;
if ~exist('bshow', 'var'), bshow=0; end;
if ~exist('brev', 'var'), brev=0; end;
if ~exist('zcal', 'var'), zcal=[]; end;
if ~exist('fldark', 'var'), fldark=[]; end;
if ~exist('ctrack', 'var'), ctrack={}; end;
if ~exist('xyscale', 'var'), xyscale=0.074; end;
if ~exist('btime', 'var'), btime=[]; end;
if ~exist('t0', 'var'), t0=0; end;
if ~exist('vframe', 'var'), vframe=[]; end;
if ~exist('vcrop', 'var'), vcrop=[]; end;
if ~exist('vresize', 'var'), vresize=[]; end;
if ~exist('nflt', 'var'), nflt=0; end;
if ~exist('tpflt', 'var'), tpflt='Linear'; end;
if ~exist('rmask', 'var'), rmask=0; end;
if ~exist('outfl', 'var'), outfl=[]; end;
if ~exist('iquad', 'var'), iquad=[]; end;
if ~exist('pivfl', 'var'), pivfl=[]; end;
if ~exist('bomg','var'), bomg=0; end;
if ~exist('zdep', 'var'), zdep=10; end;
if ~exist('npres', 'var'), npres=0; end;
if ~exist('mskrt', 'var'), mskrt=1; end;
if ~exist('mskfl', 'var'), mskfl=[]; end;
if ~exist('brel', 'var'), brel=0; end;
if ~exist('bomgz', 'var'), bomgz=0; end;
if ~exist('outflfit', 'var'), outflfit=[]; end;
if ~exist('stpvaromg', 'var'), stpvaromg=[]; end;
if ~exist('indsel', 'var'), indsel=[]; end;
end




function [datacl] = cleardup (data)
[sd, si]=sort(data(:,1));
data_srt=data(si, :);
stat_dup=(diff(data_srt(:,1))==0);
if sum(stat_dup)>0
[labeled, n]=bwlabel(stat_dup, 4);
ncurr=size(data_srt, 1);
for i=n:-1:1
inddup=find(labeled==i);
ndupi=numel(inddup);
data_srt(inddup(1):ncurr-ndupi, :)=data_srt(inddup(1)+ndupi:ncurr, :);
ncurr=ncurr-ndupi;
end
datacl=data_srt(1:ncurr, :);
else
    datacl=data;
end
end
 
