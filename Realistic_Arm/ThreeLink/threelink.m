%% Initialize dynamic model of 3-link arm
clear;
close all;

%% Robot parameters
m = [1,1,1]; m1 = m(1); m2 = m(2);m3 = m(3);
l = [1,1,1]; l1 = l(1); l2 = l(2);l3 = l(3);
c = l/2;     c1 = c(1); c2 = c(2); c3 = c(3); 
I = 0.5.*m.*l.^2; I1 = I(1); I2 = I(2);I3 = I(3);

%% Dynamic model for the two links with 4 springs
%
% system dimension
n = 3;

% state variables
x = sym('x',[2*n,1]);
q = x(1:n);
qdot = x((n+1):(2*n));

%input, which are 4 spring coefficients
u = sym('u',[2*n,1]);
k1 = u(1);
k2 = u(2);
k3 = u(3);
k4 = u(4);
k5 = u(5);
k6 = u(6);

%theta is used to compute the position of every link
theta=[q(1);q(1)+q(2);q(1)+q(2)+q(3)];
thetadot=[qdot(1);qdot(1)+qdot(2);qdot(1)+qdot(2)+qdot(3)];

%set the distance parameters for the 4 springs
d=[1,1];d1=d(1);d2=d(2); 
s=[0.4;0.4;0.3;0.4;0.3];s1=s(1);s2=s(2);s3=s(3);s4=s(4);s5=s(5);

%set the rest length of the four springs
r_rest=[0.9;0.9;0.5;0.5;0.7;0.7];
r1_rest=r_rest(1);
r2_rest=r_rest(2);
r3_rest=r_rest(3);
r4_rest=r_rest(4);
r5_rest=r_rest(5);
r6_rest=r_rest(6);

%calculate the length of two springs w.r.t q in the first link
r1=sqrt(s1^2+d1^2-2*s1*d1*cos(q(1))); 
r2=sqrt(s1^2+d2^2+2*s1*d2*cos(q(1)));
%calculate the length of two springs w.r.t q in the second link
r3=sqrt(s2^2+s3^2-2*s2*s3*sin(q(2)));
r4=sqrt(s2^2+s3^2+2*s2*s3*sin(q(2)));
%calculate the length of two springs w.r.t q in the third link
r5=sqrt(s4^2+s5^2-2*s5*s4*sin(q(3)));
r6=sqrt(s4^2+s5^2+2*s5*s4*sin(q(3)));

%%get the potential and kinetic energy for the system
pc(2,1,3) = sym('0');
pc(:,:,1) = c1*[cos(theta(1));sin(theta(1))];
pc(:,:,2) = l1*[cos(theta(1));sin(theta(1))] + ...
            c2*[cos(theta(2));sin(theta(2))];%set the Cartesian coordinates of center of mass of link
pc(:,:,3) = l1*[cos(theta(1));sin(theta(1))] + ...
            c2*[cos(theta(2));sin(theta(2))]+...
            c3*[cos(theta(3));sin(theta(3))];
%Kinetic energy
T=0;
for i = 1:n
    J = jacobian(pc(:,:,i),q); % temporary variable
    T = T + 0.5*m(i)*transpose(qdot)*( simplify(transpose(J)*J) )*qdot + ...
            0.5*I(i)*thetadot(i)^2;
    T = simplify(T);
    clear J;
end
%Potential energy
    g = 0.5; % gravity acceleration
    G=m1*g*pc(2,:,1)+m2*g*pc(2,:,2)+m3*g*pc(2,:,3);%gravity of the links
    U=.5*k1*((r1-r1_rest)^2)+.5*k2*((r2-r2_rest)^2)...
        +.5*k3*((r3-r3_rest)^2)+.5*k4*((r4-r4_rest)^2)+...
    .5*k5*((r5-r5_rest)^2)+.5*k6*((r6-r6_rest)^2);%potential energy for springs
    V=G+U;
    
 %%%get the M C phi given the lagranian expression
%fric=0;
fric=0.5*qdot;
%P=matlabFunction(fric);

%x=[q;qdot;k1;k2;k3;k4];%variables for this system
x=[q;qdot];
sym_lagrangian_solver(T,fric,x);%get the M C from the Lagrangian
%[M,C,phi_p] = Symb_lagrange_solver_varying_k(T,V_p,x);
%phi_func_p=matlabFunction(phi_p);

phi=transpose(simplify(jacobian(V,q)));
%phi_func=matlabFunction(phi);
x=[q;qdot;u];
matlabFunction(phi,  'vars',{x},'file','phi_fun');

phi_jac_x = jacobian(phi,[q;qdot]);
matlabFunction(phi_jac_x,'vars',{x},'file','phi_jac_fun_x');

phi_jac_u = jacobian(phi,u);
matlabFunction(phi_jac_u,'vars',{x},'file','phi_jac_fun_u');

%% *Simulation for targrt reaching*
dt = 1e-2; % sample time
N = 700; % total steps
t = 0:dt:((N-1)*dt);


% cost function parameters
nx=6;
nu=6;
x_init = [pi/2;0;0;0;0;0];
x_target = [-pi/5;0;0;0;0;0];
u_init = [5; 5;4;4;4;4];

Q = zeros(nx, nx, N); % initialization
R = zeros(nu, nu, N);
L = zeros(nu, N);
h = zeros(nx, N);
g = zeros(nx, N);
l_m = zeros(nu, N);
h_m = zeros(nx, N);

for k = 1:N % weighted matrix for the input u 
     R(:, :, k) = 1*eye(nu, nu);
end

Q(:, :, N) = 100000*eye(nx, nx); % final cost for the target reaching
h(:, N) = - Q(:, :, N)*x_target; 

page=3; 
x = zeros(nx, N,page);x_dev = zeros(nx, N,page); 
u = zeros(nu, N,page);u_dev = zeros(nu, N,page); 

x(:, 1) = x_init;
u(:, 1) = u_init;

A = zeros(nx, nx, N);
B = zeros(nx, nu, N);

for kk=1:page
    switch kk
        case 1 %the first linear model we just linear around the initial state

        for k=1:N
            [ A(:, :, k),B(:, :, k) ] = linearized_model( x(:, 1),u(:, 1) );
            A(:, :, k) = A(:, :, k)*dt + eye(nx);
            B(:, :, k) = B(:, :, k)*dt;    
        end
            [K, s, P, p, c] = LQR_design( A, B, Q, R, L, h, g);
    
        for k = 1:(N-1)
            u(:, k+1,kk)=K(:,:,k)*x(:, k,kk)+s(:,k);
            x(:, k+1,kk)=(A(:,:,k)+B(:,:,k)*K(:,:,k))*x(:, k,kk)+B(:,:,k)*s(:, k);
        end
        
        otherwise 
        for k=1:N   %then we keep adding derivations to the previous trajectory
            [ A(:, :, k),B(:, :, k) ] = linearized_model( x(:, k,kk-1),u(:, k,kk-1) );
            A(:, :, k) = A(:, :, k)*dt + eye(nx);
            B(:, :, k) = B(:, :, k)*dt;
            l_m(:, k) = L(:, k) + R(:, :, k)*u(:, k,kk-1);
            h_m(:, k) = h(:, k) + Q(:, :, k)*x(:, k,kk-1);
        end
            [K, s, P, p, c] = LQR_design( A, B, Q, R, l_m, h_m, g);
    
        for k = 1:(N-1)
            u_dev(:, k,kk) = K(:, :, k)*x_dev(:, k,kk) + s(:, k);
            u(:, k,kk) = u(:, k,kk-1) + u_dev(:, k,kk); % modified input
            x_dev(:, k+1,kk)=(A(:,:,k)+B(:,:,k)*K(:,:,k))*x_dev(:, k,kk)+B(:,:,k)*s(:, k);
            x(:, k+1,kk) = x(:, k+1,kk-1) + x_dev(:, k+1,kk);%modified trajectory
        end
        
       
    end
     
end
%play the movie
    for kk=1:3
      D=Threelink_vert2D_createmovie_reaching(x(:,:,kk),t,10,l,x_target);      
    end