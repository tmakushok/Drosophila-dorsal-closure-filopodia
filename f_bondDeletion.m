function [Idel] = f_bondDeletion(I, inter, RatioThres, BW, ThresLinkLen, Mask)
%% Parameters
PtNb = 15; % Number of points on the bond that are taken to do the linear fit
FillDist = 3; % The distance to the middle of the bond at which the two filling points are going to be taken
%InThresurf = 600;
OutThresurf = 400;
%%
Idel=I;     % Final watershed image will have some lines deleted
s_BW = size(BW);    % Size of the image
%% Take off double appearances of bonds between cells
% Constructing a matrix with lengths of each bond
LElemI = zeros(length(inter), 2);
for i = 1:length(inter)
    LElemI(i,:) = size(inter{i});   
end
LElemI = LElemI(:,1);
% Finding locations of double appearances
L = length(LElemI);
ToDel = [];
for i = 1:L
    ind = find(LElemI(i+1:L) == LElemI(i)) + i;
    for j = 1:length(ind)  % In most cases this loop will be executed only ones
        if inter{i} == inter{ind(j)}
            ToDel = [ToDel; ind(j)];
        end
    end
end
% Deleting the double appearances
ToDel = sort(ToDel, 'descend'); % Sort in descending order
for i = 1:length(ToDel)
    inter(ToDel(i)) = [];
end
%% Prepare the mask
Mask = ~Mask * 2;
interL = length(inter);     
Avg = zeros(interL, 1);
%% Loop for bonds deletion
for i = 1:interL    % Loop on all bonds found in the image
    Link = inter{i};
    LinkLen = size(Link,1);     % Not the length of the line, but its number of pixels
%% Get rid of the lines that touch the border of the image
    if Link(1,1) == 2 || ~isempty(find(Link(:,2) == 2)) || Link(LinkLen,1) == 2 || Link(1,1) == s_BW(1)-1 || ~isempty(find(Link(:,2) == s_BW(2)-1)) || Link(LinkLen,1) == s_BW(1)-1
        Idel(sub2ind(size(Idel), Link(2:LinkLen - 1, 1), Link(2:LinkLen - 1, 2))) = 1;  % = 1 because of inverted image
%         figure, imshow(Idel);
%         pause(0.3);
%         close
        continue
    end
%% Not analyze short lines (not to create small holes on the LE)
    if LinkLen < ThresLinkLen
       continue 
    end   
%% Find the middle coordinates in the matrix Link: 
     Middlept = Link(round(size(Link, 1)/2), 1:2);
%% Linear fit of the central piece of the bond
    % Distance between Middlept and each point of the bond
    D = sqrt((Middlept(1) - Link(:,1)) .^ 2 + (Middlept(2) - Link(:,2)) .^ 2); 
    D = [(1:size(D,1))', D];
    % Sort elements of the matrix in an ascending manner
    D = sortrows(D,2);
    % Take PtNb of the smallest elements
    LinePt = D(1:PtNb, 1);
    LinePt = [Link(LinePt, 2), Link(LinePt, 1)];   % LinePt is in the shape of [X, Y]
    % Sort points according to x position and, for the same x, according to y position
    LinePt = sortrows(LinePt);
    % If the line is vertical, then take the two points on a horizontal line
%     figure, imshow(Idel);   % To visualise the points found
%     hold on
    if LinePt(1,1) == LinePt(PtNb,1) 
        x1_Fill = Middlept(2) + FillDist;
        x2_Fill = Middlept(2) - FillDist;        
        y1_Fill = Middlept(1);
        y2_Fill = Middlept(1);
    else
        % Linear fit
        p = polyfit(LinePt(:,1), LinePt(:,2), 1);       % Result: coefficients of linear fit
        % Create the points of the line
        x1 = LinePt(1,1);                    % First X of the line
        x2 = LinePt(size(LinePt, 1),1);      % Last X of the line
        x_Line = x1 : (x2-x1)/PtNb : x2;     % All Xs of the line
        y_Line = p(1) * x_Line + p(2); 
        % Visualisation of the linear fit of the central piece of the bond        
%         plot(x_Line, y_Line, 'r*');
%         plot(Middlept(2), Middlept(1), 'o');
    %     plot([x1, x2], [p(1) * x1 + p(2), p(1) * x2 + p(2)], 'ro'); % Two ending points
    %% Finding two points for filling: they are on both sides of the middle, on perpendicular line        
        % Finding the point on the line that corresponds to the middle of the bond
        D = sqrt((Middlept(2) - x_Line) .^ 2 + (Middlept(1) - y_Line) .^ 2); 
        [a, ind] = min(D);
        x0 = x_Line(ind);
        y0 = y_Line(ind);     
        % Coeeficients of the perpendicular line k1*x + k2 = 0
        k1 = - 1/p(1);      % Relation of slopes of perpendicular lines
        k2 = y0 - k1*x0;
        % Calculating positions of the two filling points (using the angle)
        Angle = atand(k1);
        dX1 = FillDist * cosd(Angle);
        dY1 = - FillDist * sind(Angle);
        x1_Fill = x0 + dX1;
        x2_Fill = x0 - dX1;
        y1_Fill = y0 - dY1;
        y2_Fill = y0 + dY1;         
    end
    % Visualisation of the two points
%     plot([x1_Fill, x2_Fill], [y1_Fill, y2_Fill], 'g*');
%     hold off;
%% Filling from the found points and measuring surfaces of overlap       
    I1 = imfill(~I, [round(y1_Fill) round(x1_Fill)]);
    I2 = imfill(~I, [round(y2_Fill) round(x2_Fill)]);
%% Erode to get rid of other lines
    se = strel('disk',1);        
    I1 = imerode(I1, se);
    I2 = imerode(I2, se);
%     figure, imshow(I1);
%     figure, imshow(I2);
%% Find overlap with mask (In = within dorsal opening)
    Sum = I1 + Mask;
    S_In1 = length(find(Sum == 1));
    S_Out1 = length(find(Sum == 3));    
    Sum = I2 + Mask;
    S_In2 = length(find(Sum == 1));
    S_Out2 = length(find(Sum == 3));    
%% Both filled surfaces are inside
    if (S_Out1 < OutThresurf) && (S_Out2 < OutThresurf) %&& (S_Out1 > InThresurf) & (S_Out2 > InThresurf)        
        Idel(sub2ind(size(Idel), Link(2:LinkLen - 1, 1), Link(2:LinkLen - 1, 2))) = 1;  % = 1 because of inverted image
    end
%     imshow(Idel);
%     close all;
end
