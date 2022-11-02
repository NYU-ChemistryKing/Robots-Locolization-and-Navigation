function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

% Initializing
x1 = uPrev(1:3,1);
x2 = uPrev(4:6,1);
x3 = uPrev(7:9,1);
x4 = uPrev(10:12,1);
x5 = uPrev(13:15,1);

% Angvel and Noise
w = angVel;
nzero = zeros(15,1);
noise = 0.0001 * randn(15,1);
Q =  eye(15);
ng = nzero(1:3,1);
na = nzero(4:6,1);
nbg = noise(7:9,1);
nba = noise(10:12,1);

g = [-0.23; -0.195; -9.829];

% G and R functions - ZYX EulerAngles

Ginv = [
cos(x2(3))/cos(x2(2)), sin(x2(3))/cos(x2(2)), 0;
-sin(x2(3)),  cos(x2(3)), 0;
cos(x2(3))*sin(x2(2))/cos(x2(2)), sin(x2(2))*sin(x2(3))/cos(x2(2)), 1;
];

Rotation = [
    cos(x2(2))*cos(x2(3)), sin(x2(1))*sin(x2(2))*cos(x2(3))-cos(x2(1))*sin(x2(3)), cos(x2(1))*sin(x2(2))*cos(x2(3))+sin(x2(1))*sin(x2(2));
    cos(x2(2))*sin(x2(3)), sin(x2(1))*sin(x2(2))*sin(x2(3))+cos(x2(1))*cos(x2(3)), cos(x2(1))*sin(x2(2))*sin(x2(3))-sin(x2(1))*cos(x2(2));
    -sin(x2(2)), sin(x2(1))*cos(x2(2)), cos(x2(1))*cos(x2(2));
    ];

% xdot = f(x,u,n)
f1 = x3;
f2 = Ginv * (w - x4 - ng);
f3 = g + Rotation * (acc - x5 - na);
f4 = nbg;
f5 = nba;

f = [f1; f2; f3; f4; f5];

% Linearlization
At =[
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 1, 0, 0,                             0,                             0,  0,                  0,                                                0,                                               0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 1, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 1,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                    (cos(x2(3))*sin(x2(2))*(w(1) - x4(1)))/cos(x2(2))^2 + (sin(x2(2))*sin(x2(3))*(w(2) - x4(2)))/cos(x2(2))^2,                                                                            (cos(x2(3))*(w(2) - x4(2)))/cos(x2(2)) - (sin(x2(3))*(w(1) - x4(1)))/cos(x2(2)), 0, 0, 0,            -cos(x2(3))/cos(x2(2)),            -sin(x2(3))/cos(x2(2)),  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                - cos(x2(3))*(w(1) - x4(1)) - sin(x2(3))*(w(2) - x4(2)), 0, 0, 0,                      sin(x2(3)),                     -cos(x2(3)),  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,    cos(x2(3))*(w(1) - x4(1)) + sin(x2(3))*(w(2) - x4(2)) + (cos(x2(3))*sin(x2(2))^2*(w(1) - x4(1)))/cos(x2(2))^2 + (sin(x2(2))^2*sin(x2(3))*(w(2) - x4(2)))/cos(x2(2))^2,                                                          (cos(x2(3))*sin(x2(2))*(w(2) - x4(2)))/cos(x2(2)) - (sin(x2(2))*sin(x2(3))*(w(1) - x4(1)))/cos(x2(2)), 0, 0, 0, -(cos(x2(3))*sin(x2(2)))/cos(x2(2)), -(sin(x2(2))*sin(x2(3)))/cos(x2(2)), -1,                  0,                                                0,                                                0;
0, 0, 0,   (sin(x2(1))*sin(x2(3)) + cos(x2(1))*cos(x2(3))*sin(x2(2)))*(acc(2) - x5(2)) + (cos(x2(1))*sin(x2(2)) - cos(x2(3))*sin(x2(1))*sin(x2(2)))*(acc(3) - x5(3)), (cos(x2(2))*sin(x2(1)) + cos(x2(1))*cos(x2(2))*cos(x2(3)))*(acc(3) - x5(3)) - cos(x2(3))*sin(x2(2))*(acc(1) - x5(1)) + cos(x2(2))*cos(x2(3))*sin(x2(1))*(acc(2) - x5(2)), - (cos(x2(1))*cos(x2(3)) + sin(x2(1))*sin(x2(2))*sin(x2(3)))*(acc(2) - x5(2)) - cos(x2(2))*sin(x2(3))*(acc(1) - x5(1)) - cos(x2(1))*sin(x2(2))*sin(x2(3))*(acc(3) - x5(3)), 0, 0, 0,                             0,                             0,  0, -cos(x2(2))*cos(x2(3)),   cos(x2(1))*sin(x2(3)) - cos(x2(3))*sin(x2(1))*sin(x2(2)), - sin(x2(1))*sin(x2(2)) - cos(x2(1))*cos(x2(3))*sin(x2(2));
0, 0, 0, - (cos(x2(1))*cos(x2(2)) + sin(x2(1))*sin(x2(2))*sin(x2(3)))*(acc(3) - x5(3)) - (cos(x2(3))*sin(x2(1)) - cos(x2(1))*sin(x2(2))*sin(x2(3)))*(acc(2) - x5(2)), (sin(x2(1))*sin(x2(2)) + cos(x2(1))*cos(x2(2))*sin(x2(3)))*(acc(3) - x5(3)) - sin(x2(2))*sin(x2(3))*(acc(1) - x5(1)) + cos(x2(2))*sin(x2(1))*sin(x2(3))*(acc(2) - x5(2)),   cos(x2(2))*cos(x2(3))*(acc(1) - x5(1)) - (cos(x2(1))*sin(x2(3)) - cos(x2(3))*sin(x2(1))*sin(x2(2)))*(acc(2) - x5(2)) + cos(x2(1))*cos(x2(3))*sin(x2(2))*(acc(3) - x5(3)), 0, 0, 0,                             0,                             0,  0, -cos(x2(2))*sin(x2(3)), - cos(x2(1))*cos(x2(3)) - sin(x2(1))*sin(x2(2))*sin(x2(3)),   cos(x2(2))*sin(x2(1)) - cos(x2(1))*sin(x2(2))*sin(x2(3));
0, 0, 0,                                                                 cos(x2(1))*cos(x2(2))*(acc(2) - x5(2)) - cos(x2(2))*sin(x2(1))*(acc(3) - x5(3)),                                                - cos(x2(2))*(acc(1) - x5(1)) - cos(x2(1))*sin(x2(2))*(acc(3) - x5(3)) - sin(x2(1))*sin(x2(2))*(acc(2) - x5(2)),                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                                                         sin(x2(2)),                               -cos(x2(2))*sin(x2(1)),                               -cos(x2(1))*cos(x2(2));
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0;
0, 0, 0,                                                                                                                               0,                                                                                                                                        0,                                                                                                                                          0, 0, 0, 0,                             0,                             0,  0,                  0,                                                0,                                                0];

Ut = [
 
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
           -cos(x2(3))/cos(x2(2)),            -sin(x2(3))/cos(x2(2)),  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                     sin(x2(3)),                     -cos(x2(3)),  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
-(cos(x2(3))*sin(x2(2)))/cos(x2(2)), -(sin(x2(2))*sin(x2(3)))/cos(x2(2)), -1,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0, -cos(x2(2))*cos(x2(3)),   cos(x2(1))*sin(x2(3)) - cos(x2(3))*sin(x2(1))*sin(x2(2)), - sin(x2(1))*sin(x2(2)) - cos(x2(1))*cos(x2(3))*sin(x2(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0, -cos(x2(2))*sin(x2(3)), - cos(x2(1))*cos(x2(3)) - sin(x2(1))*sin(x2(2))*sin(x2(3)),   cos(x2(2))*sin(x2(1)) - cos(x2(1))*sin(x2(2))*sin(x2(3)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,           sin(x2(2)),                               -cos(x2(2))*sin(x2(1)),                               -cos(x2(1))*cos(x2(2)), 0, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                            0,                             0,  0,                  0,                                                0,                                                0, 0, 0, 0, 0, 0, 1, 0, 0, 0];


% Discretization
I = eye(15);
Ft = I + dt * At;
Vt = Ut;
Qd = dt * Q;

% Prediction
uEst = uPrev + dt * f;
covarEst = Ft * covarPrev * Ft.' + Vt * Qd * Vt.';

end


