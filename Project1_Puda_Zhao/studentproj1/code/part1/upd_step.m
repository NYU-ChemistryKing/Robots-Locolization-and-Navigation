function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively

% Observation Model
Ct = eye(6,15);
Wt = eye(6);
R =  0.01 * eye(6);
g = Ct * uEst;

% Kalman Gain
Kt = covarEst * Ct.' * pinv(Ct * covarEst * Ct.' +  Wt * R * Wt.');

% Update step
uCurr = uEst + Kt * (z_t - g) ;
covar_curr = covarEst - Kt * Ct * covarEst;

end