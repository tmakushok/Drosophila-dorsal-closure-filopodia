function [BkGd] = f_naturalImageBkGd_WithFit(InImage) 
%--------------------------------------------------------------------------
p_SmSpline = 0.09;
XOutLess = 1;
%--------------------------------------------------------------------------
[m, n] = size(InImage);
Pixels = reshape(InImage, 1, m * n);     % Making a vector out of the matrix
%% Creating a histogram of values of intensity of all pixels in an image
% NbBins = 3000;
MaxIntens = max(max(InImage));
[Nb,XOut] = hist(Pixels(find(Pixels ~= 0)), 500); % "~= 0" because croping produced pixels with intensity 0
% figure, bar(XOut, Nb); 
%% Smoothing spline applied to the histogram
% Creating an array with x values for smoothing spline dots positionning:
% dots are 10 times more frequent than the dots in the histogram analysed
x_Spl = XOut(1):((XOut(2) - XOut(1))):XOut(length(XOut));    
FitOptions = fitoptions('Method', 'SmoothingSpline', 'SmoothingParam', p_SmSpline);
FitType = fittype('smoothingspline');
cfun = fit(XOut', Nb', FitType, FitOptions);
% hold on;
% plot(x_Spl, cfun(x_Spl), '-r');
% hold off;
a = cfun(x_Spl);
[a, BkGd_index] = max(a');      
BkGd = x_Spl(BkGd_index);
     