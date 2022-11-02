function [position, orientation, R_c2w] = estimatePose(data, t)
%% CHANGE THE NAME OF THE FUNCTION TO estimatePose
% Please note that the coordinates for each corner of each AprilTag are
% defined in the world frame, as per the information provided in the
% handout. Ideally a call to the function getCorner with ids of all the
% detected AprilTags should be made. This function should return the X and
% Y coordinate of each corner, or each corner and the centre, of all the
% detected AprilTags in the image. You can implement that anyway you want
% as long as the correct output is received. A call to that function
% should made from this function.
    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    
    % Calibration
    cameraMatrix = [311.0520, 0, 201.8724; 0, 311.3885, 113.6210; 0, 0, 1];
    cameraCalibrationErrors =  [-0.04, 0.0, -0.03]; 
    Yaw = pi/4;
    % Transfrom from camera frame to imu(body) frame
    R_ci = [cos(Yaw), -sin(Yaw), 0; -sin(Yaw), -cos(Yaw), 0; 0, 0, -1;];
    T_ci = [R_ci, cameraCalibrationErrors.';0, 0, 0, 1];
    % Get pc and pw
    id = data(t).id;
    num = length(id);
    pc = zeros(2*num,4);
    pw = zeros(2*num,4);
    A = zeros(8*num,9);
    for i = 1:num
        % Get point_camera
        p1c = data(t).p1(:,i);
        p2c = data(t).p2(:,i);
        p3c = data(t).p3(:,i);
        p4c = data(t).p4(:,i);
        point_c = [p1c,p2c,p3c,p4c];
        pc(2*i-1:2*i,:) = point_c;
        % Get point_world
        point_w = getCorner(id(i));
        pw(2*i-1:2*i,:) = point_w;
        A(8*i-7:8*i,:) = getA(point_w, point_c);   
    end
    % Compute H by SVD
    H = solveH(A);
    % Solve Rotation and Translation
    B = cameraMatrix \ H;
    R1_est = B(:,1);
    R2_est = B(:,2);
    T_est = B(:,3);
    R_matrix = [R1_est, R2_est, cross(R1_est,R2_est)];
    % The second SVD to solve R & T
    [U,~,V] = svd(R_matrix);
    S1 = eye(3);
    S1(3,3) = det(U * V);

    Rotation = U * S1 * V;
    Translation = T_est / norm(R1_est);
    T_cw = [Rotation, Translation; 0,0,0,1;];
    T_wi = T_cw \ T_ci;

    Rotation = T_cw(1:3,1:3);
    %% Output Parameter Defination
    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order ZYX
    position = T_cw(1:3,4).';
    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order ZYX
    orientation = rotm2eul(Rotation,'ZYX');
    
    %R_c2w = Rotation which defines camera to world frame
    R_c2w = Rotation;
end

function H = solveH(A)
    [~,~,V] = svd(A);
    h = V(:,9);
    H = reshape(h,3,3)';
    H = H * sign(V(9,9));
end

function A = getA(P1, P2)
    P1 = P1.'; % repreing points in world frame 4*2 Matrix
    P2 = P2.'; % representing points in camera frame 4*2 Matrix
    x = P1(:, 1);
    y = P1(:, 2);
    X = P2(:, 1);
    Y = P2(:, 2);
    A = zeros(length(x(:))*2,9);
    for i = 1:length(x(:))
        a = [x(i),y(i),1];
        b = [0 0 0];
        c = [X(i);Y(i)];
        d = -c*a;
        A((i-1)*2+1:(i-1)*2+2,1:9) = [[a b;b a] d];
    end    
end