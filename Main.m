clear;      
close all
%----------------------------------
% Define folder of experiment and datasets to be analysed
ImFolder = '_InputImages/';
% Define the start frame of the movie to be analysed
Imstart = 1;
% Threshold for image segmentation (can be defined automatically based on the image histogram)
Thres = 60;
% Minimal area of objects kept on binary thresholded image
AreaThres = 20;
% Minimal number of pixels of a bond that can be deleted with 'bonddeletion'
% function (to avoid creating small holes in the LE)
ThresLienLen = 50;
% Parameters for watershed
MountainScale = 7;
PropWhite = 0.7; % for bond deletion step: proportion of white pixels on thresholded binary image on the bond
% For angle measurement
PtNb = 15;          % Better to have this number odd, to have line points symmetric around central point 
%% Filtering background images (if they were not filtered already)
f_filterBkGdImages;	
% List names of all the datasets
movielist = dir([ImFolder '*Exp_*']);
for i_mov = 1:length(movielist)         % loop on experiments 
    imagelist = dir([ImFolder movielist(i_mov).name]);    
    Maxprojlist= dir([ImFolder, movielist(i_mov).name, '/MAX_*']);
    FinTracks = cell(length(Maxprojlist), 1);            
    for i_im = Imstart:length(Maxprojlist)          % loop on images    
        close all
%% Open image
        Path = [ImFolder, movielist(i_mov).name, '/', Maxprojlist(i_im).name];
        Image = load(Path);
        Image = Image.MaxImage;
%% Illumination correction and background subtraction (for the current image)
 		Image = f_IlluminationCorrection(Image);
%% Cropping according to position of dorsal opening on the previous image
        if i_im == Imstart 
            % 1) Manual cropping for image 1 (to remove the outside of the dorsal opening)             
            figure, imshow(Image, []);
            % 'CropRect' is a four-element position vector[xmin ymin width height] 
            % that specifies the size and position of the crop rectangle
            % !!! Draw cropping rectangle at some distance from the LE !!!
            [Image, CropRect] = imcrop;
            % 2) Defining the mask for near-LE position
            % !!! Click close to the LE !!!
            Mask = f_LEdrawing(Image);       % zeros(size(Image));
%             figure, imshow(Mask, []);
            % 3) Applying Mask to Image
            Image = Mask .* Image;
        else
            % 1) Crop image according to the manual cropping performed for the
            % first image 
            Image = imcrop(Image, CropRect); 
            % 2) Mask for near-LE position is defined by previous image
            % 3) Dilate the mask a bit
            se = strel('disk', 10, 0);
            Mask = imdilate(Mask,se);            
            % 4) Applying Mask to Image
            Image = Mask .* Image;
        end   
%% Binarise the image: thresholding
        % Create binary image (using thresholding)
        BW = zeros(size(Image));
        BW(Image > Thres) = 1;
%         figure, imshow(BW, []);     
        % Take off small objects from the binary image
        BW = bwareaopen(BW, AreaThres);
%         figure, imshow(BW, []);
%         figure, imshow(Image, []);
%% Finding leading edge borders (and non-useful elements) using watershed     
		[WShed, WithAllBonds] = f_LEdetection(BW, MountainScale, PropWhite, ThresLienLen, Mask);
        figure, imshow(WithAllBonds);
        figure, imshow(WShed);
%% Finding perimeter of the leadind edge (for one image) and cropping the image (Mask)
        % Grayscale - image with only filopodia in the dorsal opening
        % Perim - matrice with LE coordinates
        % Mask - white dorsal opening 
        close all;
		[Grayscale, Perim, Mask, MaskCont] = f_perimeterDetection(WShed, Image);
        figure, imshow(Grayscale, []);
%% Detection of filopodia
		Skelet = f_Filopodiadetection(Grayscale);     % use Imagefilo as file name if subtracting amnioserosa cells
%% Linking filopodia to the LE, finding the angle between filopodia and LE		
        FilFin = f_LinkToLE(Skelet, MaskCont);        % Result: structure with fields {Real; Linked; PtLE}
%% Do visualisation: overlay only connected filopodia on top of initial image
        % Do RGB image out of out initial image
        Im = Image / max(max(Image));
        ImRGB(:,:,1) = Im; 
        ImRGB(:,:,2) = Im;
        ImRGB(:,:,3) = Im;
%         figure, imshow(ImRGB, []);
        % Find positions where connected filopodia are located
        LinesIm = zeros(size(Image));
        for i_Show = 1:length(FilFin)
            Line = round(FilFin(i_Show,1).Linked);
            LinesIm(sub2ind(size(LinesIm), Line(:,2), Line(:,1))) = 1;
        end        
        LinesIm = logical(LinesIm);
        figure, imshow(LinesIm, []);
        % Make connected filopodia in red on the RGB image
        % R:
        I = ImRGB(:,:,1);
        I(LinesIm) = 1;
        ImRGB(:,:,1) = I;
        % G:
        I = ImRGB(:,:,2);
        I(LinesIm) = 0;
        ImRGB(:,:,2) = I;
        % B:
        I = ImRGB(:,:,3);
        I(LinesIm) = 0;
        ImRGB(:,:,3) = I;
        % Show the RGB
        figure, imshow(ImRGB, []);
%% Finding angles between filopodia and LE and adding it as a new field to 'FilFin'
% !!! Carefull: Axis Oy goes down, so angles are measured in the opposite
% way from the usual.
        FilFin = f_Angle(FilFin, Perim, PtNb); 
%% Putting the results coming from different frames together
        FinTracks(i_im) = {FilFin};        
    end     % (end of loop on images)    
end % (end of loop on datasets)


