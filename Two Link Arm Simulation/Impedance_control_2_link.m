%% Impedance control of 2-link arm
% move the arm to an arbitrary position with impedance control
% must run dynamic_model_2_link.m first
% a typical test case is:
% no gravity
% friction coefficient: beta = 2

close all;

% Simulation parameters
dt = 1e-3; % sample time
N = 10000; % total steps
t = 0:dt:((N-1)*dt);
x_init = [pi/2; pi/20; 0; 0;];
x_target = [pi/2; pi/2; 0; 0];
targetCartesian = polar2Cartesian_fun(x_target(1), x_target(2));

%% Simulation Initialization

% torque
u = zeros(nu, N);

% trajectory
x = zeros(nx, N);

% initial states
x(:, 1) = x_init;

%% Run the robot arm
for k = 1:(N-1)
    % arbitrary input
    % u(:,k) = [10, 2];
    
    % input that balances gravity
%     u(:,k) = V_jac_fun(x(1,k),x(2,k));
%     
%     % input by impedance control
%     K = 8;
%     u(:,k) = u(:,k) + transpose(polar2Cartesian_jac_fun(x(1,k),x(2,k)))*K*(targetCartesian - polar2Cartesian_fun(x(1,k), x(2,k)));
    
%     % input saturation
%     u1max = 0.6;
%     u2max = 0.2;
%     if u(1,k) > u1max
%         u(1,k) = u1max;
%     elseif u(1,k) < - u1max
%         u(1,k) = - u1max;
%     end
%     if u(2,k) > u2max
%         u(2,k) = u2max;
%     elseif u(2,k) < - u2max
%         u(2,k) = - u2max;
%     end
       
    % increment of states
    dx34 = f34_fun(x(1,k), x(2,k), x(3,k), x(4,k), u(1,k), u(2,k));
    
    % calculate the state deviation of the next step
    x(1, k+1) = x(1, k) + x(3, k)*dt;
    x(2, k+1) = x(2, k) + x(4, k)*dt;
    x(3, k+1) = x(3, k) + dx34(1)*dt;
    x(4, k+1) = x(4, k) + dx34(2)*dt;
end

robot_arm_movie([l1; l2], x, x_target, 0, 10);