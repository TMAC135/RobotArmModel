clear, clc, close all;

% generate 2-link arm model synbolically
n = 2;
syms x1 x2 x3 x4
x = [x1;x2;x3;x4];
q = [x1;x2];
qdot = [x3;x4];

% arm parameters
m1 = 0.1; % mass of link 1
m2 = 0.1; % mass of link 2
l1 = 1; % length of link 1
l2 = 1; %length of link 2
r1 = l1/2; % center of mass location for link 1
r2 = l2/2; % center of mass location for link 2
Iz1 = m1*l1^2/12;
Iz2 = m2*l2^2/12;
Iz = [Iz1, 0; 
        0, Iz2];
l = [l1; l2];
m = [m1; m2];
r = [r1; r2];

% center of mass position
pc1 = [r1*cos(x1); r1*sin(x1)];
pc2 = [l1*cos(x1)+r2*cos(x1+x2); l1*sin(x1)+r2*sin(x1+x2)];
pc = [pc1, pc2];
pc1_jac = jacobian(pc1,q);
pc2_jac = jacobian(pc2,q);

% joint and tip position
joint = [l1*cos(x1); l1*sin(x1)];
tip = [l1*cos(x1)+l2*cos(x1+x2); l1*sin(x1)+l2*sin(x1+x2)];
tip_jac = jacobian(tip, q);
matlabFunction(joint, 'vars', {q}, 'file', 'joint_fun');
matlabFunction(tip, 'vars', {q}, 'file', 'tip_fun');
matlabFunction(tip_jac, 'vars', {q}, 'file', 'tip_jac_fun');

% Kinetic energy
T = 0.5*(transpose(cumsum(qdot))*Iz*cumsum(qdot) + ...
         m1*transpose(qdot)*(simplify(transpose(pc1_jac)*pc1_jac))*qdot + ...
         m2*transpose(qdot)*(simplify(transpose(pc2_jac)*pc2_jac))*qdot);
T = simplify(T);
                               
% Potential energy
g = 1;
V = m1*g*pc1(2) + m2*g*pc2(2);

V_jac = jacobian(V,q);
matlabFunction(V_jac, 'vars', {q}, 'file', 'V_jac_fun');

% Friction
beta = 1;
fric = beta*qdot;

% Generate system dynamics via Lagrangian solver
sym_lagrangian_solver(T,V,fric,x);