function [Vel] = velocityRANSAC(optV,optPos,Z,R_c2w,e)
%% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
    %% Input Parameter Description
    % optV = The optical Flow
    % optPos = Position of the features in the camera frame 
    % Z = Height of the drone
    % R_c2w = Rotation defining camera to world frame
    % e = RANSAC hyper parameter
    Flag = 1;
    if(Flag == 1)        
        num = length(optPos);
        x = optPos(1,:);
        y = optPos(2,:);
        x_dot = optV(:,1);
        y_dot = optV(:,2);
        X_A = x';
        Y_A = y';
        Z_A = Z';
        A_all = [-1./Z_A, zeros(num, 1), X_A./Z_A, X_A.*Y_A, -(1+X_A.^2), Y_A;
        zeros(num, 1), -1./Z_A, Y_A./Z_A, (1+Y_A.^2), -X_A.*Y_A, -X_A];
        ssd = @(x, y) (x-y).^2;
       
        max_inliers = zeros(num, 1);
        for t = 1:20
            ids = randi([1,num], 6, 1);
            A = Comp_A(x(ids)',y(ids)',Z(ids)');   
            V_C = A\[x_dot(ids);y_dot(ids)];
            ssd_vec = ssd([x_dot;y_dot], A_all*V_C);
            ssd_vec = ssd_vec(1:num) + ssd_vec(num+1:end);
            inliers = ssd_vec < e;
            if sum(inliers) > sum(max_inliers)
                max_inliers = inliers;
            end
        end
        A = Comp_A(x(max_inliers)',y(max_inliers)',Z(max_inliers)');
        Vel = A\[x_dot(max_inliers);y_dot(max_inliers)];
    end
end


function [A] = Comp_A(x,y,Zc)
    n = numel(x);
    A = [-1./Zc, zeros(n, 1), x./Zc, x.*y, -(1+x.^2), y;
          zeros(n, 1), -1./Zc, y./Zc, (1+y.^2), -x.*y, -x];
end