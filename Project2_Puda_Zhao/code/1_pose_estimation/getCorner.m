function res = getCorner(id)
%% CHANGE THE NAME OF THE FUNCTION TO getCorner
    %% Input Parameter Description
    % id = List of all the AprilTag ids detected in the current image(data)
    % From id to values of row and column
    row = mod(id,12);
    column = (id - row) / 12;
    gap = 0.152;
    gaperror = 0.026;
    % From row and column to coordinate of p4
    p4_y = 2 * column * gap;
    if column >= 6 
        p4_y = p4_y + 2 * gaperror;
    elseif column >= 3
        p4_y = p4_y + gaperror;
    end
    p4_x = 2 * row * gap;
    p4 = [p4_x; p4_y];
    
    % From p4 to other p points
    p3 = p4 + [0; gap];
    p1 = p4 + [gap; 0];
    p2 = p4 + [gap; gap];
    %% Output Parameter Description
    % res = List of the coordinates of the 4 corners (or 4 corners and the
    % centre) of each detected AprilTag in the image in a systematic method
    res = [p1, p2, p3, p4];
end