%% PROJECT 2 VELOCITY ESTIMATION
close all;
clear all;
clc;
addpath('../data')

%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 1;
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
estimatedV = zeros(6,length(sampledData));

%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION
% Parameters
cameraMatrix = [311.0520, 0, 201.8724; 0, 311.3885, 113.6210; 0, 0, 1];

LastF = sampledData(1).img;
LastTime = sampledData(1).t;
for n = 2:length(sampledData)
    %% Initalize Loop load images
    CurrF = sampledData(n).img;

    %% Detect good points
    points = detectFASTFeatures(LastF,'minContrast',100/255,'minQuality',30/255);
    locs = points.Location;

    %% Initalize the tracker to the last frame.
    tracker = vision.PointTracker('MaxBidirectionalError',3);
    initialize(tracker, locs, LastF);
    
    %% Find the location of the next points;
    oldpoints = locs;
    [points, isFound] = step(tracker, CurrF);

    %% Calculate velocity
    % Initialization
    CurrTime = sampledData(n).t;
    timegap = CurrTime - LastTime;
    oldpoints_C = zeros(3,length(oldpoints));
    points_C = zeros(3,length(oldpoints));
    
    % Use a for loop
    for i=1:length(oldpoints)
        oldpoints_P = [oldpoints(i,1);oldpoints(i,2);1];
        points_P = [points(i,1);points(i,2);1];
        oldpoints_C(:,i) = cameraMatrix \ oldpoints_P;
        points_C(:,i) = cameraMatrix \ points_P;
        dot_p = (points_C' - oldpoints_C') ./ timegap;
    end
    p_dot = dot_p(:,1:2);

    %% Calculate Height
    [position, orientation, R_c2w] = estimatePose(sampledData, n-1);
    Tc_w = position.';
    u = points(:,1);
    v = points(:,2);
    m = numel(u);

    C_lam = cameraMatrix \ [u';v';ones(1,m)];
    r1 = R_c2w' * Tc_w;
    r1 = r1(3);
    r2 = R_c2w' * C_lam;
    r2 = r2(3,:);
    lam = r1./r2;

    Xc = C_lam(1,:).*lam;
    Yc = C_lam(2,:).*lam;
    Zc = C_lam(3,:).*lam;
    
    %% RANSAC    
    % Write your own RANSAC implementation in the file velocityRANSAC
    result = velocityRANSAC(p_dot, oldpoints_C, Zc, R_c2w, 0.001);
    %% Thereshold outputs into a range.
    % Not necessary
    
    %% Fix the linear velocity
    % Change the frame of the computed velocity to world frame
    R_wc = R_c2w.';
    Adjoint_wc = [R_wc, zeros(3,3); zeros(3,3), R_wc];
    result = Adjoint_wc * result;

    %% ADD SOME LOW PASS FILTER CODE
    % Not neceessary but recommended 
    %estimatedV(:,n) = Vel;

    %% STORE THE COMPUTED VELOCITY IN THE VARIABLE estimatedV AS BELOW
    %estimatedV(:,n) = Vel; % Feel free to change the variable Vel to anything that you used.
    % Structure of the Vector Vel should be as follows:
    % Vel(1) = Linear Velocity in X
    % Vel(2) = Linear Velocity in Y
    % Vel(3) = Linear Velocity in Z
    % Vel(4) = Angular Velocity in X
    % Vel(5) = Angular Velocity in Y
    % Vel(6) = Angular Velocity in Z
    estimatedV(:,n) = result;
    
    LastF = sampledData(n).img;
    LastTime = CurrTime;
end

for p = 1:3
    wpass = 0.01;
    x = estimatedV(p,:);
    estimatedV(p,:) = lowpass(x, wpass);
end

for p = 4:6
    wpass = 0.15;
    x = estimatedV(p,:);
    estimatedV(p,:) = lowpass(x, wpass);
end

plotData(estimatedV, sampledData, sampledVicon, sampledTime, datasetNum)
