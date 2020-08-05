%% Finding angles between filopodia and LE and adding it as new field to 'FilFin'
function FilFin = Angle(FilFin, Perim, PtNb)
% Perim - matrix with LE coordinates
dX = 7;
for i_fil = 1:length(FilFin)    % i_fil = current filopodium number
    Real = FilFin(i_fil).Real;
    Linked = FilFin(i_fil).Linked;
    PtLE = FilFin(i_fil).PtLE;      
    % PtLE is in the shape of [X, Y] Perim is in the shape of [row col]
%% Create a LE tangent line
    % Distance between PtLE and each point of the LE perimeter
    D = sqrt((PtLE(1) - Perim(:,2)) .^ 2 + (PtLE(2) - Perim(:,1)) .^ 2);
    D = [(1:size(D,1))', D];
    % Sort elements of the matrix in an ascending manner
    D = sortrows(D,2);
    % Take PtNb of the smallest elements
    LinePt = D(1:PtNb, 1);
    LinePt = [Perim(LinePt, 2), Perim(LinePt, 1)];   % LinePt is in the shape of [X, Y]
    % Sort points according to x position
    LinePt = sortrows(LinePt);
    % Line
    p = polyfit(LinePt(:,1), LinePt(:,2), 1);       % Result: coefficients of linear fit
%% Angle measurement
    % Create two vectors from the vertices.
    % v1 = [x1 - x2, y1 - y2]
    % v2 = [x3 - x2, Y3 - y2]
    % Filopodium vector
    xV1 = FilFin(i_fil).LineEnd(1);
    yV1 = FilFin(i_fil).LineEnd(2);
    v1 = [xV1 - PtLE(1), yV1 - PtLE(2)]; 
    % LE vector  
    xV2 = PtLE(1) + dX;
    yV2 = p(1) * xV2 + p(2);
    v2 = [xV2 - PtLE(1), yV2 - PtLE(2)]; 
    % Find the angle
    theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
    % Convert it to degrees.   
    FilFin(i_fil).angle = theta * (180/pi);    
end