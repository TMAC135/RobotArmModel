syms q qdot;%m=1,l=1
m=1;
l=1;
k=1;
x=[q;qdot];
T=1/2*m*9.8*l^2*qdot^2;
V=m*9.8*l*sin(q);
fric=k*qdot;
sym_lagrangian_solver(T,V,fric,x);
%%
x_initial=[-pi/2;0];
u_initial=0;
Q=[10 0;0 9];
R=0;
N=1000;
dt=0.1;
l=1;
x_target=[pi;0];%destination point
[X,U]=mpc_control(x_initial,x_target,u_initial,Q,R,N,dt);
close all;
for T=1:N

plot([0,1*cos(X(1,T))],[0,1*sin(X(1,T))],'b',cos(x_target),sin(x_target),'b*');
axis([-l,l,-l,l]);
daspect([1,1,1]);
pause(0.1);
end
