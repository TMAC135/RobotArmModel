%%The following code is the simple code version to test the model 
%with two springs

%%
clear all;

%%the simple arm with two spring attached in the link
 syms q qdot L;
 syms g l m I r r1 r2 k d c s; 
 

%% Robot parameters
m = 1;
l = 1; 
c = l/2;  
I = 0.5.*m.*l.^2;
g=0;
k=5;
d=1;
s=0.8*c;
x_init = [pi/4; 0];
%format long;
r1_rest=sqrt(s^2+d^2);
%r1_rest=d;
%r1_rest=sqrt(s^2+d^2-2*s*d*cos(x_init(1)));

r2_rest=sqrt(s^2+d^2);
%r2_rest=sqrt(3)*d;
%r2_rest=sqrt(s^2+d^2+2*s*d*cos(x_init(1)));

r1=sqrt(s^2+d^2-2*s*d*cos(q));
%r1=sqrt(s^2+d^2-2*s*d*sin(q));
r2=sqrt(s^2+d^2+2*s*d*cos(q));
%r2=sqrt(s^2+d^2+2*s*d*sin(q));

% Cartesian coordinates of center of mass of links
pc(2,1) = sym('0');
pc(:,:) = c*[cos(q);sin(q)];
%pc(:,:) = c*[sin(q);cos(q)];

% Kinetic energy

    J = jacobian(pc(:,:),q); % temporary variable
    T = 0.5*m*transpose(qdot)*( simplify(transpose(J)*J) )*qdot + ...
            0.5*I*qdot^2;
    T = simplify(T);
%T=0.5*m*[(c)^2*(sin(q))^2*(qdot)^2+(c)^2*(cos(q))^2*(qdot)^2]+0.5*I*(q)^2;
% potential energy    
%V=m*g*c*sin(q)+k*[2*s^2+2*d^2-r1_rest*(r1+r2)];
V=m*g*c*sin(q)+.5*k*[(r1-r1_rest)^2+(r2-r2_rest)^2];

%fric=0;
beta =0.05*eye(1);
fric = beta*qdot;
 x=[q;qdot];
 %T=0.5*m*[(c)^2*(sin(q))^2*(qdot)^2+(c)^2*(cos(q))^2*(qdot)^2]+0.5*I*(q)^2;
 %V=m*g*l*sin(q)+0.5*k*[2*(l)^2+2*(d)^2+2*(r)^2-2*r*(r1+r2)];
 
[M,C,phi] = Symb_lagrange_solver(T,V,x);
%sym_lagrangian_solver(T,V,fric,x);
%M_fun(0)
%fric_fun((0,0))
Q=matlabFunction(phi);%Q is the matlabfunction of phi
P=matlabFunction(fric);%P is the matlabfunction of fric
%%here is the simulation when we add no input into the system

%% Simulation parameters
dt = 1e-2; % sample time
N = 2000; % total steps
t = 0:dt:((N-1)*dt);

nx=2;
nu=1;
% cost function parameters 
% format long;
%x_init = [pi/3; 0];
u = zeros(nu, N);
x = zeros(nx, N);
x(:, 1) = x_init;
q=x_init(1,1);
qdot=x_init(2,1);
%GET F and G for every iteration
for k = 1:(N-1)
        F=[0 ,1;0,0];
        %G=[0;M\(-Q(q))];
        G=[0;M\(-Q(q))+M\(-P(qdot))];
        %G=[0;4*(-Q(q))/3+4*(-G(qdot))/3];
        %G=[0;M\phi_fun(x(:,k))];
        %G=[0;M_fun(q)\phi_fun(q)+M_fun(q)\fric_fun(x(:,k))];
        F=F*dt+eye(nx);
        G=G*dt;
 
        x(:, k+1) = F*x(:, k)+G;
        q=x(1,k+1);
        qdot=x(2,k+1);
        
end
figure;
plot(t,x(1,:),'r');
legend('angle');

hold on;
plot(t,x(2,:),'b');
hold off;
%legend('angular velocity');
D=threelink_vert2D_createmovie(x,t,10,l);