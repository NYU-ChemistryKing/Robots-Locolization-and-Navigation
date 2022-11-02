clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);

Z = sampledVicon(1:6,:);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.1*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0; %last time step in real time
pos = proj2Data.position;
pose = proj2Data.angle;
z_t = [pos.'; pose.'];
for i = 1:length(sampledTime)
    %% Fill in the FOR LOOP
    currTime = sampledTime(i); 
    %Predction Step
    [covarEst,uEst] = pred_step(uPrev, covarPrev, sampledData(i).omg, sampledData(i).acc, currTime - prevTime);
    % Update Step
    [uCurr,covar_curr] = upd_step(z_t(:,i), covarEst, uEst);
    % Save current step 
    savedStates(:,i) = uCurr;
    uPrev = uCurr;
    covarPrev = covar_curr;
    % Update time
    prevTime = currTime;
end

plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);