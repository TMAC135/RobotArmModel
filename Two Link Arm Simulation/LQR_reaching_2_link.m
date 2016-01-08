%% LQR controled two-link planar arm reaching
% This program uses iterative LQR algorithm to find a locally optimal 
% solution of a reaching task
% The robustness of the algorithm highly depends on the friction of the
% system. Up to now (02/26/2015), the best way found to increase robustness
% is by introducing friction.
% must run dynamic_model_2_link.m first
% typical test case:
% m1 = 1, l1 = 1, m2 = 0.5, l2 = 0.5, beta = 1
% robust reaching could be performed to every target

close all;

% Simulation parameters
dt = 1e-3; % sample time
N = 10000; % total steps
t = 0:dt:((N-1)*dt);

% cost function parameters
x_init = [0; 0; 0; 0];
x_target = [-pi; 9*pi/10; 0; 0];
Q = zeros(nx, nx, N); % initialize
R = zeros(nu, nu, N);
l = zeros(nu, N);
h = zeros(nx, N);
g = zeros(nx, N);
for k = 1:N % intermidiate cost
     R(:, :, k) = 0.001*eye(nu, nu);
end
Q(:, :, N) = 1e5*eye(nx, nx); % final cost
h(:, N) = - Q(:, :, N)*x_target;

% modified cost function parameters for linearized model
l_m = zeros(nu, N);
h_m = zeros(nx, N);

%% Simulation Initialization

% torque
u = zeros(nu, N);
for ii = 1:N
    u(:,ii) = V_jac_fun(0,0);
end
u_dev = zeros(nu, N); % deviation

% trajectory in each iteration
x = zeros(nx, N);
x_dev = zeros(nx, N); % deviation

% initial states
x(:, 1) = x_init;

% Linearized model in each iteration
A = zeros(nx, nx, N);
B = zeros(nx, nu, N);

% nonlinear model
D = zeros(nu, nu);
C = zeros(nu, nx);

% total cost
cost_LQR = 0;

% calculate the partial derivatives used in the linearized model
% partial_derivatives;

%% iterations

for kk = 1:10 % total number of iterations
    % modify model and cost
    for k = 1:N
        % linearized system model around the previous trajectory
        par_ders = f34_jac_fun( x(1,k), x(2,k), x(3,k), x(4,k), u(1,k), u(2,k) );     
        A(:, :, k) = [0 0 1 0;
                      0 0 0 1;
                      par_ders(1:2, 1:4)]*dt + eye(nx, nx);
        B(:, :, k) = [0 0;
                      0 0;
                      par_ders(1:2, 5:6)]*dt;
        % modified cost function parameters
        l_m(:, k) = l(:, k) + R(:, :, k)*u(:, k);
        h_m(:, k) = h(:, k) + Q(:, :, k)*x(:, k);
    end
    % design LQR controler
    [K, s, P, p, c] = LQR_design( A, B, Q, R, l_m, h_m, g);
    
%     % total cost
%     if kk == 1 
%         cost_LQR = transpose(x_init)*P(:,:,1)*x_init + ...
%                transpose(p(:,1))*x_init + c(1);
%     else
%         cost_LQR = cost_LQR + c(1);
%     end
    
    % run the robot arm
    for k = 1:(N-1)
        % calculate the input deviation by the control law
        u_dev(:, k) = K(:, :, k)*x_dev(:, k) + s(:, k);
        u(:, k) = u(:, k) + u_dev(:, k); % modified input
                
        % drive the arm with input
        dx34 = f34_fun(x(1,k), x(2,k), x(3,k), x(4,k), u(1,k), u(2,k));
        
        % calculate the state deviation of the next step
        x_dev(1, k+1) = x(1, k) + x(3, k)*dt - x(1, k+1);
        x_dev(2, k+1) = x(2, k) + x(4, k)*dt - x(2, k+1);
        x_dev(3, k+1) = x(3, k) + dx34(1)*dt - x(3, k+1);
        x_dev(4, k+1) = x(4, k) + dx34(2)*dt - x(4, k+1);
        
        % calculate the state of the next step
        x(:, k+1) = x(:, k+1) + x_dev(:, k+1);
    end
    
    if kk > 0
        robot_arm_movie([l1; l2], x, x_target, 0, 1); % the last inuput is the rough time of the movie
    end
        
%     % calculate cost by definition
%     cost_def = 0;
%     for k = 1:N
%         cost_def = transpose(x(:,k))*Q(:,:,k)*x(:,k) + transpose(u(:,k))*R(:,:,k)*u(:,k) + ...
%                      2*transpose(l(:,k))*u(:,k) + 2*transpose(h(:,k))*x(:,k) + cost_def;
%     end
    
%     % plot results
%     figure,
%     plot(t,x(1,:),t,x(2,:),t,x(3,:),'-.',t,x(4,:),'-.',t,u(1,:),'--',t,u(2,:),'--');
%     xlabel('Time (seconds)'), ylabel('x1,x2: rad    x3,x4: rad/s'), grid;
%     legend('Theta1','Theta2','Theta1 dot','Theta2 dot','Torque1','Torque2'),
%     title({'2-link arm reaching'; ...
%            ['cost by def = ', num2str(cost_def), ' cost by LQR = ', num2str(cost_LQR)]; ...
%            ['initial state = ', mat2str(x_init), ' target = ', mat2str(x_target)]});
    
end