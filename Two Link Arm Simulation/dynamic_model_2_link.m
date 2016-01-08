%% Create Dynamic Model of 2-link Arm for Multiple Purposes
% different potential energies can be specified
% friction is included
% most of the useful functions are written to files

clear, clc;
close all;

%% Robot parameters
m1 = 1; % mass of link 1
m2 = 0.5; % mass of link 2
l1 = 1; % length of link 1
l2 = 0.5; %length of link 2
r1 = l1/2; % center of mass location for link 1
r2 = l2/2; % center of mass location for link 2

% Define alpha, beta and delta from masses and lengths (for simulation)
Iz1 = m1*l1^2/12;
Iz2 = m2*l2^2/12;
a = Iz1 + Iz2 + m1*(r1^2) + m2*(l1^2 + r2^2);
b = m2*l1*r2;
d = Iz2 + m2*(r2^2);

% system dimension
nx = 4;
% input dimension
nu = 2;

%% Dynamic model with symbols

syms x1 x2 x3 x4 tau1 tau2
syms x tau
x = [x1; x2; x3; x4];
tau = [tau1; tau2];

gravity = 1e1; % gravity acceleration
beta = 2; % friction coefficient

% Potential Energy
G = r1*sin(x1)*m1*gravity + ( l1*sin(x1) + r2*sin(x1+x2) )*m2*gravity; 
JointSpring = 1e1*heaviside(x2 - 9*pi/10)*(x2 - 9*pi/10)^2 + 1e1*heaviside(-x2-0.2)*(-x2-0.2)^2;
V = G; % V coule be the combination of the above different potential energies

V_jac_syms = simplify(jacobian(V, [x1;x2]));
V_jac_handle = matlabFunction(V_jac_syms, 'vars', [x1, x2, x3, x4], 'file', 'V_jac_fun.m');

% Kinetic Energy
% because the 2-link arm case is relatively simple, the D and C matrices
% are written by hand in here. In the future, this will be replaced by the
% more general function that generate the dynamic model from Lagrangian
D_syms = [ a + 2*b*cos(x2), d+b*cos(x2);
          d + b*cos(x2)  ,      d      ];
C_syms = [ - 2*b*sin(x2)*x4, -b*sin(x2)*x4;
             b*sin(x2)    ,       0       ];
D_handle = matlabFunction(D_syms, 'vars', x, 'file', 'D_fun.m');
C_handle = matlabFunction(C_syms, 'vars', x, 'file', 'C_fun.m');

% dynamic equations
f34_syms = simplify(D_syms\( - C_syms*[x3; x4] + tau - transpose(V_jac_syms) - beta*[x3; x4]));
f34_handle = matlabFunction(f34_syms, 'vars', [x; tau], 'file', 'f34_fun.m');
f34_jac_syms = simplify(jacobian(f34_syms, [x; tau]));
f34_jac_handle = matlabFunction(f34_jac_syms, 'vars', [x; tau], 'file', 'f34_jac_fun.m');

%% Coordinates transformation
polar2Cartesian_syms = [l1*cos(x1) + l2*cos(x1+x2); l1*sin(x1) + l2*sin(x1+x2)];
polar2Cartesian_handle = matlabFunction(polar2Cartesian_syms, 'vars', [x1, x2], 'file', 'polar2Cartesian_fun.m');
polar2Cartesian_jac_syms = jacobian(polar2Cartesian_syms);
polar2Cartesian_jac_handle = matlabFunction(polar2Cartesian_jac_syms, 'vars', [x1, x2], 'file', 'polar2Cartesian_jac_fun.m');