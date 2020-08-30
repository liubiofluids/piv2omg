function [zdepth]=r2zdepth(R)
load('r2zdepthdat.mat');
zdepth=interp1(vR, vmsz, R);
end