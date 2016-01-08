model_2_link_arm;
x_initial = [-1; -1; 0; 0];
x_target = [-0.9;0; 0; 0];
u_initial=[0;0];
Q=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]*1000;
R=[1 0;0 1];
l=[1;1];
N=1000;
dt=0.1;

[X,U]=mpc_control(x_initial,x_target,u_initial,Q,R,N,dt);

%%
robot_arm_movie_2_link( l, X(1:2,:), x_target, 1000);