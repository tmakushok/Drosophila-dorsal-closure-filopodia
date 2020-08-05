function [Skelet] = f_Filopodiadetection(InitImage)
%% Parameters
%!!!--!!! 'Rigidity' parameter for smoothing spline
p_OneMTFit = 0.001;  % for smoothing before calculating total bundle length per cell
p_AllMTFit = 0.05; %0.005; for maxima
%!!!--!!! Threshold for BkGd subtraction on LoG-filtered images
BkGdThres = - 0.03;       % 0.8;  
%!!!--!!! Min skeleton length to be considered as a bundle
MinSkeletLen = 8;
%%
BundlesLengths = [];        % To store length of each bundle with its number and number of the cell
TotalMTLengthAllCells = [];
TotalMTLengthAllCells_PerCellLengthUnit = [];
FigNb = 1;
UntreatedIm = InitImage;
InitFigNb = figure, imshow(InitImage, []);
%% Using Laplacian of Gaussian for smoothing filopodia preserving their tips 
h = fspecial('log', 7, 4);  
LapGaus = imfilter(InitImage, h, 'replicate'); 
LapGaus(find(LapGaus > BkGdThres)) = 0;
% figure, imshow(LapGaus, []);
BinaryIm = zeros(size(LapGaus));
BinaryIm(find(LapGaus < 0)) = 1;   
% figure, imshow(BinaryIm, []);   
%% Skeletonization
Skelet = bwmorph(BinaryIm, 'thin', Inf);
s = size(Skelet);
%% Adding lines at the bottom and right
Skelet = [Skelet, zeros(s(1), 1); zeros(1, s(2) + 1)];
% figure, imshow(Skelet);
%% Chosing long fragments of the skeleton
%Converting binary image to a label matrix   
Labels = bwlabel(Skelet);
Stats = regionprops(Labels, 'Area'); 
% Taking off very small labelled fragments of the skeleton
Skelet = ismember(Labels, find([Stats.Area] > MinSkeletLen));    
% figure, imshow(Skelet);
% Creating a black frame around the image (to avoid having filopodia
% detection on the sides of the image)
Skelet(1,:) = 0;
Skelet(s(1),:) = 0;
Skelet(:,1) = 0;
Skelet(:,s(2)) = 0;
%% Overlay skeleton on top of initial cropped image
[SkelLin, SkelCol] = find(Skelet);
figure, imshow(UntreatedIm, []); hold on;
for i = 1:length(SkelLin)    
    line(SkelCol(i), SkelLin(i), 'Color', [.8 0 0], 'Marker', '.'); 
end
%% Finding all filopodia ends used later for 3D tracking
% Calculating an image where at a previously white on 'Skelet' image pixel 
% there is the sum of the pixel and of all surrounding pixels
WhiteSkelPos = find(Skelet);        % Linear indexes of all white pixels
s = size(Skelet);
ImMTEnds = zeros(s);
for i = 1:length(WhiteSkelPos)    
    [i_Lin, j_Lin] = ind2sub(s, WhiteSkelPos(i));
    ImMTEnds(i_Lin, j_Lin) = Skelet(i_Lin - 1, j_Lin - 1) + Skelet(i_Lin - 1, j_Lin) + ...
        Skelet(i_Lin - 1, j_Lin + 1) + Skelet(i_Lin, j_Lin - 1) + Skelet(i_Lin, j_Lin) + ...
        Skelet(i_Lin, j_Lin + 1) + Skelet(i_Lin + 1, j_Lin - 1) + Skelet(i_Lin + 1, j_Lin) + ...
        Skelet(i_Lin + 1, j_Lin + 1);
end
%% Filopodia tips found
% figure, imshow(ImMTEnds, []);  
[MTTipsI MTTipsJ] = find(ImMTEnds == 2);






