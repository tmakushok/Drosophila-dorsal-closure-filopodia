%PERIMETER DETECTION
function [Grayscale, PerimCoords, Mask, MaskCont] = f_perimeterDetection(S, Grayscale)
%% Center of the image
s = size(S);
CP = [floor(s(2)/2), floor(s(1)/2)];    % [x, y]
%% Find positions where values are black and make them white
Sinvert = ~S;
%% Dilate image
se = strel('disk',1,0);
Sdil = imdilate(Sinvert,se);
% figure, imshow(Sdil, []);
%% Fill in the central object
Mask = imfill(Sdil, [CP(2) CP(1)], 4);
% figure, imshow(Mask, []);
%% Erode to get rid of the lines around
se = strel('disk',2,0);
Mask = imerode(Mask,se);
% figure, imshow(Mask);
%% Dilate ( 1 pixels)
se = strel('disk',1,0);
Mask = imdilate(Mask,se);
% figure, imshow(Mask);
%% Multiply initial grayscale image x Mask
Grayscale = double(Mask).*Grayscale;
% figure, imshow(Grayscale, []);
%% Find perimeter
MaskCont = bwperim(Mask,4);
% figure, imshow(MaskCont);
%% Find perimeter coordinates
[row, col] = find(MaskCont);
PerimCoords = [row, col];








%% Backup
% 

%S = watershed image
% S = imread(WatIm);
% image(S);
% imagesc(S); 

% % To find the dorsal opening
% Labels = bwlabel(Sperim);
% % To find center of all objects
% Stats = regionprops(Labels, 'Centroid');
% 
% %nb= number of objects detected
% nb = size(Stats);
% % preallocation
% D = zeros(nb(1),1); 
% 
% 
% for i = 1:nb(1)    %loop on all objects detected
%     %D = Distance between each center and the choosen central point of the object of interest
%     X = Stats(i).Centroid(1);
%     Y = Stats(i).Centroid(2);
%     
%     D(i) = sqrt((X -CP(1)) .^ 2 + (Y-CP(2)) .^ 2);
% end
% 
% % Find the smallest value in the matrice D
% [a,Obj] = min(D);
% 
% % To get positions of the perimeter of the object
% [row, col] = find(Labels == Obj);
% 
% % Create a new matrice full of "0" and add the object in
% MaskCont = zeros(size(Labels));
% MaskCont(find(Labels == Obj))=1;
% figure, imshow(MaskCont, []);

% % Fill in the object
% Mask = imfill(MaskCont,4,'holes');
% figure, imshow(Mask, []);




% %Get initial grayscale image
% %double: changes the type of matrice (bit)
% Path = ['.\Initimages\' expname '\' moviename '\' moviename '_1.tif'];
% Grayscale = double(imread(Path));
% figure, imshow(Grayscale, []);
