%% *Simulation for LQR targrt reaching*
dt = 1e-2; % sample time
N = 1000; % total steps
t = 0:dt:((N-1)*dt);

% simulation parameters
n=2;
nx=4;
nu=4;
u_init = [0;0;0;0];
x_init = [pi/2;0;0;0];
x_target = [pi/3-0.4;pi/6;0;0];
qddot = zeros(n,N);

% cost function parameters
Q = zeros(nx, nx, N); % initialize
R = zeros(nu, nu, N);
l = zeros(nu, N);
h = zeros(nx, N);
g = zeros(nx, N);
for k = 1:N % intermidiate cost
     R(:, :, k) = 0.1*eye(nu, nu);
end
Q(:, :, N) = 1e5*eye(nx, nx); % final cost
h(:, N) = - Q(:, :, N)*x_target;

% Linearized model around the initial condition
A = zeros(nx, nx, N);
B = zeros(nx, nu, N);
for k = 1:N
    % linearized system model around the previous trajectory
    [ A(:,:,k), B(:,:,k) ] = linearized_model_real_arm(x_init,u_init);
    A(:, :, k) = A(:, :, k)*dt + eye(nx);
    B(:, :, k) = B(:, :, k)*dt;
%     g(:,k)=dt*[ [zeros(nx/2),eye(nx/2)]*x_init; qddot_fun_real_arm(x_init,u_init) ];
end

% design lQR controller
[K, s, P, p, c] = LQR_design( A, B, Q, R, l, h, g);

% start simulation
u = zeros(nu, N);
x = zeros(nx, N);
x(:,1) = x_init;

% run simulation
for k = 1:(N-1)
    
    u(:,k) = K(:, :, k)*x(:, k) + s(:, k);
    
    qddot(:,k) = qddot_fun_real_arm(x(:,k),u(:,k));
    x(:,k+1) = x(:,k) + dt*[x(n+1:2*n,k);qddot(:,k)];
            
end
close all;
real_arm_movie_2_link(length,x,x_target,50,r_rest,d,ss);