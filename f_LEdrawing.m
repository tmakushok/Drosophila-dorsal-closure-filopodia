%% Drawing the leading edge by hand using "ginput"
% ginput gathers an unlimited number of points until you press the Enter key
function [MaskLE] = f_LEdrawing(Image)
figure, imshow(Image, []);
% Select points for Mask *make sure it is closed (fish shape) so it gives 1 big hole*
[x,y] = ginput;
MaskLE = zeros(size(Image));
for i_Line = 2:length(x)
    % Getting coordinates of points along the current line
    [cx,cy,c] = improfile(Image, [x(i_Line-1), x(i_Line)], [y(i_Line-1), y(i_Line)]);
    for i_pt = 1:length(cx) 
        MaskLE(int16(cy(i_pt)), int16(cx(i_pt))) = 1;
    end
end
% Dilate the contour before filling
MaskLE = imdilate(MaskLE, strel('disk', 3, 0));
% Fill
MaskLE = imfill(MaskLE, 'holes');
% Erode the mask to get rid of the 'tails' due to the 1st and last clicks
MaskLE = imerode(MaskLE, strel('disk', 4, 0));

