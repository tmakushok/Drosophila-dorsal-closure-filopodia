%% Linking filopodia to the leading edge
function [Res] = f_LinkToLE(Skelet, MaskCont)
%% Parameters
% Thresholds for lengths (in pixels) of the skeleton objects 
FilLenMin = 10;
FilLenMax = 250;
% Maximal distance (in pixels) from a line end to LE at which linkage is still allowed
MinDist2LE = 3;
MaxDist2LE = 27;
% Distance (in pixels) from a point on the continuation of the line to LE
% to see if the line does the link to LE 
%%
Res = struct('Real', {}, 'Linked', {}, 'PtLE', {});
%% Take off LE line, make all filopodia unlinked to LE 
se = strel('disk', 2, 0);
MCD = imdilate(MaskCont,se);  
% figure, imshow(MCD, []);
L = bwlabel(Skelet, 8);
s = size(L);
Lsepar = L(1:s(1)-1, 1:s(2)-1);   % There is one row and one column too much in L in comparison to MaskCont
Lsepar(find(MCD)) = 0;
% figure, imshow(Lsepar, [-20 50]);
%% Taking off lines that are longer or shorter than normal filopodia
Lsepar = bwlabel(Lsepar, 8);
StatsArea = regionprops(Lsepar, 'Area');
ind = find(([StatsArea.Area] > FilLenMin) .* ([StatsArea.Area] < FilLenMax));
Lsepar = ismember(Lsepar, ind);
% figure, imshow(Lsepar, []);
Lsepar = bwlabel(Lsepar, 8);
Stats = regionprops(Lsepar, 'PixelList');
% Prepare visualisation for filopodia checking
Separ = zeros(size(Lsepar));
Separ(find(Lsepar)) = 1;
figure, imshow(Separ + 2*MaskCont, []), hold on;
%% Linking each line to the LE contour
[LE(:,2), LE(:,1)] = find(MaskCont);    % X and Y of LE contour
for i_S = 1:length(Stats)               % loop on all lines
    % List of pixels belonging to the current line. First column = x(columns), second = y (rows)
    LinePx = Stats(i_S).PixelList;    
    s_LP = size(LinePx);  
%% Linear approximation of the line
    p = polyfit(LinePx(:,1), LinePx(:,2), 1);
    LinApp = p(1) * LinePx(:,1) + p(2);
    % Going from max line_plus length to max X_plus
    XPlus = ceil(sqrt(MaxDist2LE ^2 / (1 + p(1) ^2)));  
    if abs(p(1)) > 1
        a = 1 : 1/abs(p(1)) : XPlus;   % To have one point per pixel instead of big distances 
    else
        a = 1:XPlus; 
    end
    % For each end of the line find distances to each point of the LE
    Dist2LE1 = sqrt((LinePx(1,1) - LE(:,1)).^2 + (LinePx(1,2) - LE(:,2)).^2);
    Dist2LE2 = sqrt((LinePx(s_LP(1), 1) - LE(:,1)).^2 + (LinePx(s_LP(1), 2) - LE(:,2)).^2);
    % For each end of the line find minimal distance to the LE
    [MinDist1, IndMin1] = min(Dist2LE1);
    [MinDist2, IndMin2] = min(Dist2LE2);
%% Preparing X for the plolonged line
    coeff = sign(p(1));
    if MinDist1 < MinDist2      % It means first line end is closer to LE                     
        b = ones(1, length(a)) * double(int16(LinePx(1,1)));
        if LinApp(1,1) > LinApp(s_LP(1),1)
            LinePlus = (b + coeff * a)';
        else
            LinePlus = (b - coeff * a)';  % int16(LinApp(1,1)) * ones(1,XPlus)
        end
    else                        % It means the other end of the line is closer to LE       
        b = ones(1, length(a)) * double(int16(LinePx(s_LP(1),1)));
        if LinApp(s_LP(1),1) > LinApp(1,1)
            LinePlus = (b + coeff * a)';
        else
            LinePlus = (b - coeff * a)';  % int16(LinApp(1,1)) * ones(1,XPlus)
        end       
    end                            
%% Prologation of the line in the direction of the line from first line end on     
    LinePlus(:, 2) = p(1) * LinePlus(:, 1) + p(2);        
%% Calculating distances between each plus_line point and LE
    for i_Pt = 1:length(LinePlus(:,1)) 
        [Dist(i_Pt, 1), LEpos(i_Pt, 1)] = min(sqrt((LinePlus(i_Pt,1) - LE(:,1)).^2 + (LinePlus(i_Pt,2) - LE(:,2)).^2));
        line(LinePlus(i_Pt,1), LinePlus(i_Pt,2), 'Marker', '*', 'MarkerEdgeColor', 'r');
    end
%% Finding interseption point between LE and continuation of the line
    [DistToLE_Fin, pos] = min(Dist);
    % Checking if the line actually does cross the LE
    if DistToLE_Fin > MinDist2LE
        continue
    end
    ThePt = LE(LEpos(pos), :);
    % Visualisation
    line(ThePt(1), ThePt(2), 'Marker', 'o');  
%% Constructing linked to LE version of the filopodium (starting with the tip, ending with LE point)
    % Organize 'LinePx' so that it starts with the tip
    if MinDist1 < MinDist2          % It means first line end is closer to LE 
        LinePx = flipud(LinePx);
    end
    FilLinked = [LinePx; LinePlus(1:pos, :)];    
%% Accumulate the results
    ind = length(Res) + 1;
    Res(ind, 1).Real = LinePx;          % Filopodium coordinates before any treatment
    Res(ind, 1).Linked = FilLinked;     % Starting with the tip, ending with LE point
    Res(ind, 1).PtLE = ThePt;           % Intersection of filopodium with LE
    Res(ind, 1).LineEnd = LinePlus(1,:);           % Point away of the LE on linear fit of filopodium (for angle measurement)
    % Cleaning up
    Dist = [];
    LEpos = [];
end
