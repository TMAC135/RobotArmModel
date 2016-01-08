%% Impedance control of 2-link arm
% move the arm to an arbitrary position with impedance control
% must run dynamic_model_2_link.m first
% a typical test case is:
% no gravity
% friction coefficient: beta = 2

clc, close all;

% Simulation parameters
dt = 1e-3; % sample time
N = 10000; % total steps
t = 0:dt:((N-1)*dt);
x_init = [-pi/2; 0; 0; 0; 0; 0;];
x_target = [pi/2; 0; pi/2; 0; 0; 0];
targetCartesian = polar2Cartesian_fun(x_target);

%% Simulation Initialization

% torque
u = zeros(n, N);

% trajectory
x = zeros(2*n, N);

% initial states
x(:, 1) = x_init;

%% Run the robot arm
for k = 1:(N-1)
%     % arbitrary input
%     u(:,k) = [0;0;0];
    
    % input that balances gravity
    u(:,k) = phi_fun(x(:,k));
    
    % input by impedance control
    K = 8;
    u(:,k) = u(:,k) + ... 
             transpose(polar2Cartesian_jac_fun(x(:,k)))*K* ...
             (targetCartesian - polar2Cartesian_fun(x(:,k)));
         % add damping
       
    % angular acceleration
    qddot = qddot_fun(x(:,k), u(:,k));
    
    % calculate the state deviation of the next step
    x(:,k+1) = x(:,k) + dt*[ [zeros(n,n),eye(n,n)]*x(:,k); qddot ];
end

robot_arm_movie_3_link(l, x, targetCartesian, 10);