% Clear variables
clear; 

% CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
datasetNum = 9;   
[sampledData, sampledVicon, sampledTime] = init(datasetNum);

%all the measurements that you need for the update
Z = sampledVicon(1:6,:);   

% Set initial condition
% Copy the Vicon Initial state
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1));
% Covariance constant
covarPrev = eye(15); 
%Just for saving state his.
savedStates = zeros(15, length(sampledTime)); 
%last time step in real time
prevTime = 0; 

%write your code here calling the pred_step.m and upd_step.m functions
for i = 1:length(sampledTime)
    
    currTime = sampledTime(i); 
%Predction Step
    [covarEst,uEst] = pred_step(uPrev, covarPrev, sampledData(i).omg, sampledData(i).acc, currTime - prevTime);
% Update Step
    [uCurr,covar_curr] = upd_step(Z(:,i), covarEst, uEst);
% Save current step 
    savedStates(:,i) = uCurr;
    uPrev = uCurr;
    covarPrev = covar_curr;
    prevTime = currTime;
end

plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);