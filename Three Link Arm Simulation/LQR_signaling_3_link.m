%% LQR controled three-link planar arm signaling
% This program uses iterative LQR algorithm to find two different
% trajectories that best differentiate between each othe whil limiting the
% reach cost to a certain level
% The robustness of the algorithm highly depends on the friction of the
% system. The best way found to increase robustness
% is by introducing friction.
% must run dynamic_model_3_link.m first

clc, close all;

% Simulation parameters
dt = 1e-3; % sample time
N = 5000; % total steps
t = 0:dt:((N-1)*dt);

% initial conditions and targets
x_init_p = [-pi/2; pi; -pi/20; 0; 0; 0];
x_init_n = x_init_p;
x_init = [x_init_p; x_init_n];
x_target_p = [pi/20; -pi/20; pi/20; 0; 0; 0];
x_target_n = - x_target_p;
x_target = [x_target_p; x_target_n];
targetCartesian_p = polar2Cartesian_fun(x_target_p);
targetCartesian_n = polar2Cartesian_fun(x_target_n);
targetCartesian = [targetCartesian_p;targetCartesian_n];

% initialize cost functions
W = 1e5; % weights on reach error
Q = zeros(4*n, 4*n, N); % initialize
R = zeros(2*n, 2*n, N);
l = zeros(2*n, N);
h = zeros(4*n, N);
g = zeros(4*n, N);
for k = 1:N % intermidiate cost
     R(:, :, k) = 0.001*eye(2*n);
end
% if penalize reach error in polar
Q(:, :, N) = W*eye(4*n); % final cost
h(:, N) = - Q(:, :, N)*x_target;
% % if penalize reach error in Cartesian
% Q(:, :, N) = W*diag([0 0 0 1 1 1 0 0 0 1 1 1]); % penalise end point speed
% h(:, N) = - Q(:, :, N)*x_target;

% modified cost function parameters for linearized model
l_m = zeros(2*n, N);
h_m = zeros(4*n, N);
Q_m = zeros(4*n, 4*n, N);

%% Simulation Initialization

% trajectory initialized as initial position
x = zeros(4*n, N);
for ii = 1:N
    x(:,ii) = x_init;
end
x_dev = zeros(4*n, N); % deviation

% torque initialized as such balancing the potential energy at initial
% position 
u = zeros(2*n, N);
for ii = 1:N
    u(:,ii) = [ phi_fun(x(1:n,ii)); phi_fun(x((n+1):2*n,ii))];
end
u_dev = zeros(2*n, N); % deviation

% Linearized model in each iteration
A = zeros(4*n, 4*n, N);
B = zeros(4*n, 2*n, N);

% total cost
cost_LQR = 0;

%% iterations
lambda = 1e3;
for kk = 1:12 % total number of iterations
    
    if kk == 2
        % if penalize reach error in Cartesian
        Q(:, :, N) = W*diag([0 0 0 1 1 1 0 0 0 1 1 1]); % penalise end point speed
        h(:, N) = - Q(:, :, N)*x_target;
    end
    
    % modify model and cost
    for k = 1:N
        % linearized system model around the previous trajectory
        [ A(      1:2*n,       1:2*n, k), B(      1:2*n,     1:n,   k)] = linearized_model(x(1:2*n,k),u(1:n,k));
        [ A((2*n+1):4*n, (2*n+1):4*n, k), B((2*n+1):4*n, (n+1):2*n, k)] = linearized_model(x((2*n+1):4*n,k),u((n+1):2*n,k));
        A(:, :, k) = A(:, :, k)*dt + eye(4*n);
        B(:, :, k) = B(:, :, k)*dt;
        
        % modified cost function parameters
        l_m(:, k) = l(:, k) + R(:, :, k)*u(:, k);
        h_m(:, k) = h(:, k) + Q(:, :, k)*x(:, k) - signaling_jac_fun(x(:,k),[l1;l2;l3])/2/lambda;
        Q_m(:,:,k) = Q(:,:,k) - signaling_jac2_fun(x(:,k),[l1;l2;l3])/2/lambda;
    end
    % if penalise reach error in Cartesian
    if kk > 1
        h_m(:, N) = h_m(:, N) + W*reach_error2_jac_fun(x(:,N),targetCartesian);
        Q_m(:,:,N) = Q_m(:,:,N) + W*reach_error2_jac2_fun(x(:,N),targetCartesian);
    end
    
    % design LQR controler
    [K, s, P, p, c] = LQR_design( A, B, Q_m, R, l_m, h_m, g);
    
    % run the robot arm
    for k = 1:(N-1)
        % calculate the input deviation by the control law
        u_dev(:, k) = K(:, :, k)*x_dev(:, k) + s(:, k);
        u(:, k) = u(:, k) + u_dev(:, k); % modified input
        
        for i = 0:1
            % drive the arm with input
            % angular acceleration
            qddot = qddot_fun(x((1+2*n*i):(2*n+2*n*i),k), u((1+n*i):(n+n*i),k));

            % calculate the state deviation of the next step
            x_dev((1+2*n*i):(2*n+2*n*i),k+1) = ...
                x((1+2*n*i):(2*n+2*n*i),k) - ...
                x((1+2*n*i):(2*n+2*n*i),k+1) + ...
                dt*[ [zeros(n,n),eye(n,n)]*x((1+2*n*i):(2*n+2*n*i),k); qddot ];
        end
        
        % calculate the state of the next step
        x(:, k+1) = x(:, k+1) + x_dev(:, k+1);
    end
    
    % plot arm movement
    if kk > 0
        robot_arm_movie_3_link([l1, l2, l3], x, x_target, 1);
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