%% *Simulation for free move*
dt = 1e-2; % sample time
N = 1000; % total steps
t = 0:dt:((N-1)*dt);

% simulation parameters
n=2;
nx=4;
nu=4;

% initial and ending condition
x_init = [pi/2;0;0;0];
x_target = [pi/20;pi/20;0;0];
qddot = zeros(n,N);

%use inverse kinematic to specify the constant input
% choosing the magnitude of stiffness is very tricky
u_init = [3;0;3;0];
u_init(2) = K2_fun(x_target,u_init);
u_init(4) = K4_fun(x_target,u_init);

% start simulation
u = zeros(nu, N);
x = zeros(nx, N);
x(:,1) = x_init;

%GET F and G for every iteration
for k = 1:(N-1)

    u(:,k) = u_init;
    qddot(:,k) = qddot_fun_real_arm(x(:,k),u(:,k));
    x(:,k+1) = x(:,k) + dt*[x(n+1:2*n,k);qddot(:,k)];
            
end
close all;
real_arm_movie_2_link(length,x,x_target,100,r_rest,d,ss);
save('x_tra','x');
save('u_tra','u')