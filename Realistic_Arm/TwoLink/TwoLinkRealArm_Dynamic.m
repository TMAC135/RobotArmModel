%%   The following code generates the Dynamics of the TwoLinkRealArm 

%% Initialize dynamic model of 2-link real arm
clear, clc;
close all;

%% Robot parameters
mass = [1;1];       mass1 = mass(1);        mass2 = mass(2);
length = [1;1];     length1 = length(1);    length2 = length(2);
clength = length/2; clength1 = clength(1);  clength2 = clength(2); 
rot_inertia = mass.*length.^2/12; Inertial1 = rot_inertia(1); Inertial2 = rot_inertia(2);

%% Dynamic model for the two links with 4 springs
%
% system dimension
n = 2;

% state variables
x = sym('x',[2*n,1]); x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
q = x(1:n);
qdot = x((n+1):(2*n));

%input, which are 4 spring stiffness coefficients
u = sym('u',[2*n,1]);
k1 = u(1);
k2 = u(2);
k3 = u(3);
k4 = u(4);

%theta is used to compute the position of each link
theta=[q(1);q(1)+q(2)];
thetadot=[qdot(1);qdot(1)+qdot(2)];

%set the distance parameters for the 4 springs
d1=1;
d2=1*d1;
d=[d1,d2];
ss1=0.4*d1;
ss2=0.4*d1;
ss3=0.3*d1;
ss=[ss1;ss2;ss3]; %we use ss to avoid the conflict of the s in LQR design

%set the rest length of the four springs
r1_rest=1;
r2_rest=r1_rest*1;
r3_rest=r1_rest*0.4;
r4_rest=r1_rest*0.4;
r_rest=[r1_rest;r2_rest;r3_rest;r4_rest];

%calculate the length of two springs in terms of q driving the first link
r1=sqrt(ss1^2+d1^2-2*ss1*d1*cos(q(1))); 
r2=sqrt(ss1^2+d2^2+2*ss1*d2*cos(q(1)));
%calculate the length of two springs in terms of q driving the second link
r3=sqrt(ss2^2+ss3^2-2*ss2*ss3*sin(q(2)));
r4=sqrt(ss2^2+ss3^2+2*ss2*ss3*sin(q(2)));
%%get the potential and kinetic energy for the system
pc(2,1,2) = sym('0');
pc(:,:,1) = clength1*[cos(theta(1));sin(theta(1))];
pc(:,:,2) = length1*[cos(theta(1));sin(theta(1))] + ...
            clength2*[cos(theta(2));sin(theta(2))];%set the Cartesian coordinates of center of mass of link

% joint and tip position
joint = [length1*cos(x1); length1*sin(x1)];
tip = [length1*cos(x1)+length2*cos(x1+x2); length1*sin(x1)+length2*sin(x1+x2)];
tip_jac = jacobian(tip, q);
matlabFunction(joint, 'vars', {q}, 'file', 'joint_fun');
matlabFunction(tip, 'vars', {q}, 'file', 'tip_fun');
        
%Kinetic energy
T=0;
for i = 1:n
    J = jacobian(pc(:,:,i),q); % temporary variable
    T = T + 0.5*mass(i)*transpose(qdot)*( simplify(transpose(J)*J) )*qdot + ...
            0.5*rot_inertia(i)*thetadot(i)^2;
    T = simplify(T);
    clear J;
end

%Potential energy
    g = 0; % gravity acceleration
    G=mass1*g*pc(2,:,1)+mass2*g*pc(2,:,2);%gravity of the links
    U=.5*k1*((r1-r1_rest)^2)+.5*k2*((r2-r2_rest)^2)...
        +.5*k3*((r3-r3_rest)^2)+.5*k4*((r4-r4_rest)^2);%potential energy for springs
    V=G+U;

%friction
fric=1*qdot;

%genereate the M,C and phi function for future use
x=[q;qdot];
sym_lagrangian_solver_real_arm(T,V,fric,x,u);%get the M C from the Lagrangian