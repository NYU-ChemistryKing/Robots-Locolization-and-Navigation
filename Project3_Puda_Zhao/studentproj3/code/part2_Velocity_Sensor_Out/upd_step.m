function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the state
    %uEst - estimated mean of the state

    %% Measurement Model
    % Parameters
    VV = 1e-4 * eye(3);
    alpha = 0.0001;
    beta = 2;
    k = 1;
    Yaw = pi/4;
    n = length(uEst);
    %Tb2c_b =  [-0.04, 0.0, -0.03];
    Tb2c_b = [0.04 * cos(Yaw); -0.04 * sin(Yaw); -0.03];
    skew_r = @(x)[0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
    v_c = z_t(1:3,1);
    w_c = z_t(4:6,1);

    % Transfrom from imu(body) frame to camera frame 
    R = @(q)eul2rotm(flipud(q).');
    R_cb = [cos(Yaw), -sin(Yaw), 0; -sin(Yaw), -cos(Yaw), 0; 0, 0, -1;];
    R_bc = R_cb.';

    % Measurement Model - function g
    g = @(x2,x3)( R_cb * R(x2)^-1 * x3 - R_cb * skew_r(Tb2c_b) * R_bc * w_c );

    %% Step1: Compute Sigma Points - Here I did't augment x
    % Decide λ
    lamda = alpha^2*(n+k)-n;

    % Get x_aug, u_aug, P_aug
    x_aug = uEst; % 18 * 1
    P_aug = covarEst; % 18 * 18

    % Go on
    S = sqrtm(P_aug); % S 18*18
    S = sqrt(lamda+n)*S;
    W = [S, -S]; % W 18*36
    X_aug = x_aug*ones(1,2*n) + W;

    %% Step2: Propagate Sigma Points
    Z = zeros(3,2*n);
    Z0 = g(uEst(4:6,1), uEst(7:9,1));
    for i=1:2*n
        Z(:,i) = g(X_aug(4:6,i), X_aug(7:9,i));
    end

    %% Step3: Compute the predicted mean and covariance
    Wm0 = lamda/(n+lamda);
    Wc0 = lamda/(n+lamda) + (1-alpha^2+beta);
    Wm = 1/2/(n+lamda);
    Wc = 1/2/(n+lamda);

    % Calculate uEst and covarEst  
    % Update Mean
    zEst = Wm*sum(Z,2) + Wm0*Z0; 
    
    % Update Covariance
    xdiff0 = uEst - uEst;
    zdiff0 = Z0 - zEst;

    Ct = Wc0 * (xdiff0 * zdiff0.');
    St = Wc0 * (zdiff0 * zdiff0.') + VV;
    
    for i=1:2*n
        xdiff = X_aug(1:15,i)-uEst; % diff 15*1
        zdiff = Z(:,i)-zEst;
        Ct = Ct + Wc * (xdiff * zdiff.'); % 15*3
        St = St + Wc * (zdiff * zdiff.');
    end
    
    % Calculate Kalman Gain
    Kt = Ct / St;
    
    % Update uCurr and covarCurr
    uCurr = uEst + Kt * (v_c - zEst);
    covar_curr = covarEst - Kt*St*Kt.';
end
