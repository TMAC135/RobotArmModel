%%
clear all;

%%the simple arm with two spring attached in the link
 syms q qdot;
 syms g l m I r r1 r2 k1 k2 d c s; %k1 and k2 are spring coefficient for 
                                   %for the two springs
 

%% arm parameters
m = 1;
l = 1; 
c = l/2;  
I = 0.5.*m.*l.^2;
g=0;
d=1;
s=0.8*c;
x_init = [pi/2; 0];
%set the rest length of the two spring
%r1_rest=sqrt(s^2+d^2);
r1_rest=0.8;
%r2_rest=sqrt(s^2+d^2);
r2_rest=0.8;
%calculate the length of two springs w.r.t q
r1=sqrt(s^2+d^2-2*s*d*cos(q));
r2=sqrt(s^2+d^2+2*s*d*cos(q));

%%get the potential and kinetic energy for the system
pc(2,1) = sym('0');
pc(:,:) = c*[cos(q);sin(q)];%set the Cartesian coordinates of center of mass of link
%Kinetic energy

    J = jacobian(pc(:,:),q); % temporary variable
    T = 0.5*m*transpose(qdot)*( simplify(transpose(J)*J) )*qdot + ...
            0.5*I*qdot^2;
    T = simplify(T);
%Potential energy
    G=m*g*pc(2,1);%gravity of the links
    U=.5*k1*[(r1-r1_rest)^2]+.5*k2*[(r2-r2_rest)^2];%potential energy for springs
    V=G+U;

 %%%get the M C phi given the lagranian expression
 
fric=0.5*qdot;
P=matlabFunction(fric);
x=[q;qdot;k1;k2];%state variable for the system,
[M,C,phi] = Symb_lagrange_solver_varying_k(T,V,x);
Q=matlabFunction(phi);


%% *Simulation for targrt reaching*
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
k1=10;k2=5;

for k = 1:(N-1)
        F=[0 ,1;0,0];
        %G=[0;M\(-Q(k1,k2,q)+M\(-P(qdot)))];
        G=[0;M\(-Q(k1,k2,q))];
        F=F*dt+eye(nx);
        G=G*dt;
 
        x(:, k+1) = F*x(:, k)+G;
        q=x(1,k+1);
        qdot=x(2,k+1);
        k1=5/N+k1;
        k2=-5/N+k2;
end

D=threelink_vert2D_createmovie(x,t,10,l);

