%% LQR controled three-link planar arm reaching
% This program uses iterative LQR algorithm to find a locally optimal 
% solution of a reaching task
% The robustness of the algorithm highly depends on the friction of the
% system. The best way found to increase robustness
% is by introducing friction.
% must run dynamic_model_3_link.m first
% robust reaching could be performed to every target

clc, close all;

% Simulation parameters
dt = 1e-2; % sample time
N = 1000; % total steps
t = 0:dt:((N-1)*dt);

% cost function parameters
x_init = [-pi/2; pi; -pi/20; 0; 0; 0];
x_target = [pi/20; pi/20; pi/20; 0; 0; 0];
targetCartesian = polar2Cartesian_fun(x_target);

Q = zeros(2*n, 2*n, N); % initialize
R = zeros(n, n, N);
l = zeros(n, N);
h = zeros(2*n, N);
g = zeros(2*n, N);
for k = 1:N % intermidiate cost
     R(:, :, k) = 0.001*eye(n, n);
end
% if penalize reach error in polar
Q(:, :, N) = 1e5*eye(2*n, 2*n); % final cost
h(:, N) = - Q(:, :, N)*x_target;
% if penalize reach error in Cartesian
% Q(:, :, N) = W*diag([0 0 0 1 1 1]); % penalise end point speed
% h(:, N) = - Q(:, :, N)*x_target;

% modified cost function parameters for linearized model
l_m = zeros(n, N);
h_m = zeros(2*n, N);

%% Simulation Initialization

% trajectory in each iteration
x = zeros(2*n, N);
for ii = 1:N
    x(:,ii) = x_init;
end
x_dev = zeros(2*n, N); % deviation

% torque
u = zeros(n, N);
for ii = 1:N
    u(:,ii) = phi_fun(x(:,ii));
end
u_dev = zeros(n, N); % deviation

% initial states
x(:, 1) = x_init;

% Linearized model in each iteration
A = zeros(2*n, 2*n, N);
B = zeros(2*n, n, N);

% total cost
cost_LQR = 0;

% calculate the partial derivatives used in the linearized model
% partial_derivatives;

%% iterations

for kk = 1:12 % total number of iterations
    % modify model and cost
    for k = 1:N
        % linearized system model around the previous trajectory
        [ A(:, :, k), B(:, :, k)] = linearized_model(x(:,k),u(:,k));
        A(:, :, k) = A(:, :, k)*dt + eye(2*n);
        B(:, :, k) = B(:, :, k)*dt;
        
        % modified cost function parameters
        l_m(:, k) = l(:, k) + R(:, :, k)*u(:, k);
        h_m(:, k) = h(:, k) + Q(:, :, k)*x(:, k);
    end
    
%     % if penalize reach error in Cartesian
%     h_m(:, N) = h(:, N) + Q(:,:,N)*x(:,N) + ...
%                 W*transpose(reach_error_jac_fun(x(:,N),targetCartesian))/2;
    
    % design LQR controler
    [K, s, P, p, c] = LQR_design( A, B, Q, R, l_m, h_m, g);
    
    % run the robot arm
    for k = 1:(N-1)
        % calculate the input deviation by the control law
        u_dev(:, k) = K(:, :, k)*x_dev(:, k) + s(:, k);
        u(:, k) = u(:, k) + u_dev(:, k); % modified input
                
        % drive the arm with input
        % angular acceleration
        qddot = qddot_fun(x(:,k), u(:,k));

        % calculate the state deviation of the next step
        x_dev(:,k+1) = x(:,k) - x(:,k+1) + ...
                       dt*[ [zeros(n,n),eye(n,n)]*x(:,k); qddot ];
        
        % calculate the state of the next step
        x(:, k+1) = x(:, k+1) + x_dev(:, k+1);
    end
    
    % plot arm movement
    if kk > 0
        robot_arm_movie_3_link([l1, l2, l3], x, x_target, 2);
    end
    
%     % total cost by LQR
%     if kk == 1 
%         cost_LQR = transpose(x_init)*P(:,:,1)*x_init + ...
%                transpose(p(:,1))*x_init + c(1);
%     else
%         cost_LQR = cost_LQR + c(1);
%     end
% 
%     % calculate cost by definition
%     cost_def = 0;
%     for k = 1:N
%         cost_def = transpose(x(:,k))*Q(:,:,k)*x(:,k) + transpose(u(:,k))*R(:,:,k)*u(:,k) + ...
%                      2*transpose(l(:,k))*u(:,k) + 2*transpose(h(:,k))*x(:,k) + cost_def;
%     end
end