%% LQR controled two-link planar robot arm
% Reaching with signaling
% balance the cost of reaching and amount of signaling through lambda
% if lambda is too small, which means signaling is much more important than
% reaching, then the movement of arms won't converge
% must run dynamic_model_2_link.m first
% typical test case:
% l1 = l2 = 1, m1 = m2 = 1, beta = 2
% ability to converge depends on dt, N and possibly other simulation 
% parameters

close all;

%% Parameters
% Simulation parameters
dt = 1e-3; % sample time
N = 1e3; % total steps
t = 0:dt:((N-1)*dt);

% cost function parameters
W = 1e5; % cost weight on reach error
x_init = [-pi/2; pi; 0; 0];
x_target_p = [pi/10; pi/10; 0.0; 0.0];
x_target_n = [-pi/10; - pi/10; 0.0; 0.0];
x_target = [x_target_p; x_target_n];

M = [1 0 0 0 -1 0 0 0;
     0 1 0 0 0 -1 0 0];

%% iterations

% The intention of introducing index ii is to simulate different lambda in
% one run of the program. 

for ii = 1:4
    lambda = 10^(2*ii)/100;
    
    % cost function initialization
    Q = zeros(2*nx, 2*nx, N); % initialize
    R = zeros(2*nu, 2*nu, N);
    l = zeros(2*nu, N);
    h = zeros(2*nx, N);
    g = zeros(2*nx, N);
    for k = 1:N % intermidiate cost
        R(:, :, k) = 0.1*eye(2*nu, 2*nu);
    end
 
    % modified cost function parameters for linearized model
    l_m = zeros(2*nu, N);
    h_m = zeros(2*nx, N);
    
    % initialize states with initial condition
    x = zeros(2*nx, N);
    x(:, 1) = [x_init; x_init];
    u = zeros(2*nu, N);
     
%     % initialize states with previously stored trajectory
%     % current stored trajectory is the solution to lambda = 1000
%     x = importdata('x_tra_init.mat');
%     u = importdata('u_tra_init.mat');
    
    u_dev = zeros(2*nu, N); % deviation
    x_dev = zeros(2*nx, N); % deviation
    
    % Linearized model in each iteration
    A = zeros(2*nx, 2*nx, N);
    B = zeros(2*nx, 2*nu, N);
    
    % nonlinear model
    D = zeros(nu, nu);
    C = zeros(nu, nu);
    
    % total cost
    cost_LQR = 0;
    
    for kk = 1:12 % total number of iterations for one value of lambda
        % modify model and cost
        for k = 1:N
            for i = 0:1
                % linearized system model around the previous trajectory
                par_ders = f34_jac_fun( x(1+nx*i,k), x(2+nx*i,k), x(3+nx*i,k), x(4+nx*i,k), u(1+nu*i,k), u(2+nu*i,k) );
                A((1+nx*i):(nx+nx*i), (1+nx*i):(nx+nx*i), k) = [0 0 1 0;
                                                                0 0 0 1;
                                                                par_ders(1:2, 1:4)]*dt + eye(nx, nx);
                B((1+nx*i):(nx+nx*i), (1+nu*i):(nu+nu*i), k) = [0 0;
                                                                0 0;
                                                                par_ders(1:2, 5:6)]*dt;
            end
            % modified cost function parameters
            J = M*[ polar2Cartesian_jac_fun(x(1,k), x(2,k)), zeros(2, 6);
                    zeros(2, 8);
                    zeros(2, 4), polar2Cartesian_jac_fun(x(5,k), x(6,k)), zeros(2, 2);
                    zeros(2, 8)                                                        ];
            Q(:, :, k) = - transpose(J)*J/lambda;
            h_m(:, k) = h(:, k) - transpose(J)*(polar2Cartesian_fun(x(1,k), x(2,k)) - polar2Cartesian_fun(x(5,k), x(6,k)))/lambda;
            l_m(:, k) = l(:, k) + R(:, :, k)*u(:, k);
        end
        % modified reach error cost
        Q(:, :, N) = W*eye(2*nx, 2*nx) + Q(:, :, N);
        h(:, N) = - W*eye(2*nx, 2*nx)*x_target;
        h_m(:, N) = h(:, N) + Q(:, :, N)*x(:,N) - transpose(J)*(polar2Cartesian_fun(x(1,N), x(2,N)) - polar2Cartesian_fun(x(5,N), x(6,N)))/lambda;
        
        % design LQR controler
        [K, s, P, p, c] = LQR_design( A, B, Q, R, l_m, h_m, g);
        
        % run simulation
        for k = 1:(N-1)
            % calculate the input deviation by the control law
            u_dev(:, k) = K(:, :, k)*x_dev(:, k) + s(:, k);
            u(:, k) = u(:, k) + u_dev(:, k); % modified input
            
%             % input saturation
%             u1max = 0.6;
%             u2max = 0.2;
%             if u(1,k) > u1max
%                 u(1,k) = u1max;
%             elseif u(1,k) < - u1max
%                 u(1,k) = - u1max;
%             end
%             if u(2,k) > u2max
%                 u(2,k) = u2max;
%             elseif u(2,k) < - u2max
%                 u(2,k) = - u2max;
%             end
%             if u(3,k) > u1max
%                 u(3,k) = u1max;
%             elseif u(3,k) < - u1max
%                 u(3,k) = - u1max;
%             end
%             if u(4,k) > u2max
%                 u(4,k) = u2max;
%             elseif u(4,k) < - u2max
%                 u(4,k) = - u2max;
%             end
            
            for i = 0:1
                dx34 = f34_fun(x(1+i*nx,k), x(2+i*nx,k), x(3+i*nx,k), x(4+i*nx,k), u(1+i*nu,k), u(2+i*nu,k));
%                 dx34 = f34_fun(x(1+i*nx,k), x(2+i*nx,k), x(3+i*nx,k), x(4+i*nx,k), 0, 0); % no input, for testing
                
                % calculate the state deviation of the next step
                x_dev(1+i*nx, k+1) = x(1+i*nx, k) + x(3+i*nx, k)*dt - x(1+i*nx, k+1);
                x_dev(2+i*nx, k+1) = x(2+i*nx, k) + x(4+i*nx, k)*dt - x(2+i*nx, k+1);
                x_dev(3+i*nx, k+1) = x(3+i*nx, k) + dx34(1)*dt - x(3+i*nx, k+1);
                x_dev(4+i*nx, k+1) = x(4+i*nx, k) + dx34(2)*dt - x(4+i*nx, k+1);
            end
            % calculate the state of the next step
            x(:, k+1) = x(:, k+1) + x_dev(:, k+1);
        end
        
        % plot arm movement
        % kk > 10 means that we do not show the result of the first 10
        % iterations. Set kk = 1:12 so the last two iterations are shown
        % and we can tell whether the solutions have converged
        % set kk > 0 to plot all iterations
        if kk > 10
            robot_arm_movie([l1; l2], x, x_target, lambda, 1);
        end
                
        %     % total cost
        %     if kk == 1
        %         cost_LQR = transpose([x_init; - x_init])*P(:,:,1)*[x_init; - x_init] + ...
        %                transpose(p(:,1))*[x_init; - x_init] + c(1);
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
        %
        %     % plot results
        %     figure,
        %     plot(t,x(1,:),t,x(2,:),t,x(3,:),'-.',t,x(4,:),'-.',t,u(1,:),'--',t,u(2,:),'--');
        %     xlabel('Time (seconds)'), ylabel('x1,x2: rad    x3,x4: rad/s'), grid;
        %     legend('Theta1','Theta2','Theta1 dot','Theta2 dot','Torque1','Torque2'),
        %     title({'2-link arm reaching'; ...
        %            ['cost by def = ', num2str(cost_def), ' cost by LQR = ', num2str(cost_LQR)]; ...
        %            ['initial state = ', mat2str(x_init), ' target = ', mat2str(x_target)]});
        
    end
end