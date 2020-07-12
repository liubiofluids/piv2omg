function [imgbl] = im2blur (img0, gsrad)
fltrad=5*gsrad;
imgsz=size(img0);
PSF=fspecial('gaussian', fltrad, gsrad);
img1=[flipud(img0); img0; flipud(img0)];
img1=[fliplr(img1), img1, fliplr(img1)];
imgbl=imfilter(img1, PSF);
imgbl=imgbl(imgsz(1)+1:2*imgsz(1), imgsz(2)+1:2*imgsz(2));
end
