%% This function is inspired by the code by Jérôme Solon
function [S, BW] = f_LEdetection(BW2, MountainScale, PropWhite, ThresLienLen, Mask)
%% Watershed
%---------------------------------------------------------------------
D = bwdist(BW2);
fim = -D;
fim2 = imhmin(fim, MountainScale);
S2 = watershed(fim2);
%% Contour detection
inter = f_simplecontour(S2);
%% Bond deletion
BW = ones(size(S2));
BW(find(S2 == 0)) = 0;
% figure, imshow(BW);
[S]= f_bondDeletion(BW, inter, PropWhite, BW2, ThresLienLen, Mask);





