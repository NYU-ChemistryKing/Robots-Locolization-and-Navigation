function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 

    % Parameters
    alpha = 0.25;
    beta = 2;
    k = 3;
    g = [-0.23; -0.2; -9.829];
    QQ = diag([ones(1,6)*1e-2, ones(1,6)*1e-2]); % process noise cov (12*12) [ng;na;nbg;nba]
    noise = 1e-1 * randn(12,1);
    n = length(uPrev) + length(noise);

    % R and G 
    % Z-Y-X Euler Angle with q = [roll_, pitch_y, yaw_x]
    %R = @(q)eul2rotm(q.','XYZ');
    R = @(q)eul2rotm(flipud(q).');
    G = @(q)[
        cos(q(1))*cos(q(2)), -sin(q(1)), 0;
        cos(q(1))*sin(q(2)), cos(q(1)), 0;
        -sin(q(1)), 0, 1;
        ];

    % Process Model
    % get function f
    f = @(x,u,n)[
        x(7:9,1);
        G(x(4:6,1)) \ (u(1:3,1)-x(10:12,1)-n(1:3,1));
        g + R(x(4:6,1))*(u(4:6,1)-x(13:15,1)-n(4:6,1));
        n(7:9,1);
        n(10:12,1);
        ];

    % Step1: Compute Sigma Points
    % Decide λ
    lamda = alpha^2*(n+k)-n;

    % Get x_aug, u_aug, P_aug
    x_aug = [uPrev; noise];
    u_aug = [angVel; acc; zeros(n-length(angVel)-length(acc),1)];
    P_aug = [covarPrev, zeros(15,12); zeros(12,15), QQ];

    % Go on
    S = chol(P_aug); % S 27*27
    S = sqrt(lamda+n)*S;
    W = [S, -S]; % W 27*54
    X_aug = x_aug*ones(1,2*n) + W;

    % Step2: Propagate Sigma Points
    Y = zeros(15,2*n);
    Y0 = f(x_aug(1:15,1), u_aug, x_aug(16:27,1))*dt + x_aug(1:15,1);
    for i=1:2*n
        Y(:,i) = f(X_aug(1:15,i), u_aug, X_aug(16:27,i))*dt + X_aug(1:15,i);
    end

    % Step3: Compute the predicted mean and covariance
    Wm0 = lamda/(n+lamda);
    Wc0 = lamda/(n+lamda) + (1-alpha^2+beta);
    Wm = 1/2/(n+lamda);
    Wc = 1/2/(n+lamda);
    
    % Calculate uEst and covarEst  
    uEst = Wm*sum(Y,2) + Wm0*Y0; 

    diff0 = Y0-uEst;
    covarEst = Wc0 * (diff0 * diff0.');

    for i=1:2*n
        diff = Y(:,i)-uEst; % diff 27*1
        covarEst = covarEst + Wc * (diff .* diff.');
    end
end
