%% The task of the program is to smooth strongly the two images needed for background correction
function [] = f_filterBkGdImages()
% Paths
ImageBkGdFile = 'WithLazer.tif';
ImageCCDBkGdFile = 'WOLazer.tif';
BkGdImage_Output = 'Filtered_WithLazer.mat';
CCDBkGdImage_Output = 'Filtered_WOLazer.mat';
%%
CCDBkGdImage = double(imread(ImageCCDBkGdFile));        % CCD camera problem image
BkGdImage = double(imread(ImageBkGdFile)); 
%% Strong filtering of both correction images
% First a bit of median filtering to get rid of some of the noise
BkGdImage = medfilt2(BkGdImage, [3 3]);
CCDBkGdImage = medfilt2(CCDBkGdImage, [3 3]);
% Then a strong average filter
h = fspecial('average', 15);      % Creating the kernel for averaging filter
BkGdImage = imfilter(BkGdImage, h, 'replicate');
CCDBkGdImage = imfilter(CCDBkGdImage, h, 'replicate');
%% Subtraction of CCD camera problem image (filtered) from the illumination correction image   
BkGdImage = BkGdImage - CCDBkGdImage;
% imshow(BkGdImage, []);
%% Normalising intensity levels of imhomogeneous lazer illumination image to 1
% Saved now because needed for the 'cropImage' function
BkGdImage = BkGdImage / max(max(BkGdImage));
%% Output 
save(BkGdImage_Output, 'BkGdImage');  
%% Cropping of CCD problem image in the same way as analysed images are going to be cropped
% CCDBkGdImage = cropImage(BkGdImage, CCDBkGdImage);  % 'BkGdImage' is used as template for cropping
%% Visualisation
% figure, imshow(BkGdImage, []);
% figure, imshow(CCDBkGdImage, [185, 193]); 
%% Output
save(CCDBkGdImage_Output, 'CCDBkGdImage');  

