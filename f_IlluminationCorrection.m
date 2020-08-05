%% The task of the program is to correct the input maximum projection image for 
%% CCD camera shading problem and for inhomogeneous lazer illumination
%% (using correction images preferably collected on the same day)
function [InitImage] = f_IlluminationCorrection(InitImage)
ImageBkGdFile = 'Filtered_WithLazer.mat';
ImageCCDBkGdFile = 'Filtered_WOLazer.mat';
%% Loading correction images
CCDBkGdImage = load(ImageCCDBkGdFile);        % CCD camera problem image
CCDBkGdImage = CCDBkGdImage.CCDBkGdImage;
BkGdImage = load(ImageBkGdFile);       
BkGdImage = BkGdImage.BkGdImage;
%% Subtraction of already prepared CCD camera problem image from the current image   
InitImage = InitImage - CCDBkGdImage; 
%% Finding the background value
BkGdValue = f_naturalImageBkGd_WithFit(InitImage);      
%% Filtering (a bit) of the image 
h = fspecial('gaussian');   % The default value for hsize is [3 3]; the default value for sigma is 0.5. 
InitImage = imfilter(InitImage, h);    
%% Subtracting the background           
InitImage = InitImage - BkGdValue;
InitImage(InitImage < 0) = 0;  
%     figure, imshow(InitImage, []);        
%% Divide by the background image
InitImage = InitImage ./ BkGdImage;                    
%     figure, imshow(InitImage, []);                                                                                         

