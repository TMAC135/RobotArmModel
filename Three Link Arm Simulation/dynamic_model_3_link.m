%% Initialize dynamic model of 3-link arm
clear, clc;
close all;

%% Robot parameters
m = [1,0.8,0.3]; m1 = m(1); m2 = m(2); m3 = m(3);
l = [1,0.8,0.3]; l1 = l(1); l2 = l(2); l3 = l(3);
c = l/2;     c1 = c(1); c2 = c(2); c3 = c(3); 
I = 0.5.*m.*l.^2; I1 = I(1); I2 = I(2); I3 = I(3);

%% Dynamics

% system dimension
n = 3;

% state variables
x = sym('x',[2*n,1]);
q = x(1:n);
qdot = x((n+1):(2*n));
theta = cumsum(q); % angle of each link w.r.t. global coordinate
thetadot = cumsum(qdot); % angular velocities of links

% Cartesian coordinates of center of mass of links
pc = sym(zeros(2,3));
pc(:,1) = c1*[cos(theta(1));sin(theta(1))];
pc(:,2) = l1*[cos(theta(1));sin(theta(1))] + ...
          c2*[cos(theta(2));sin(theta(2))];
pc(:,3) = l1*[cos(theta(1));sin(theta(1))] + ...
          l2*[cos(theta(2));sin(theta(2))] + ...
          c3*[cos(theta(3));sin(theta(3))];

% Kinetic energy
T = 0;
for i = 1:n
    J = jacobian(pc(:,i),q); % temporary variable
    T = T + 0.5*m(i)*transpose(qdot)*( simplify(transpose(J)*J) )*qdot + ...
            0.5*I(i)*thetadot(i)^2;
    T = simplify(T);
    clear J;
end

% Potential energy
grvt = 1; % gravity acceleration
V = m1*grvt*pc(2,1) + m2*grvt*pc(2,2) + m3*grvt*pc(2,3);
V = simplify(V);
V = 0;

% friction
beta = 4*eye(n,n);
fric = beta*qdot;

% Solve for the dynamics symbolicly and generate functions
Symb_lagrange_solver(T,V,fric,x);

%% Coordinates transformation
polar2Cartesian = l1*[cos(theta(1));sin(theta(1))] + ...
                  l2*[cos(theta(2));sin(theta(2))] + ...
                  l3*[cos(theta(3));sin(theta(3))];
matlabFunction(polar2Cartesian, 'vars', {x}, 'file', 'polar2Cartesian_fun.m');
polar2Cartesian_jac = jacobian(polar2Cartesian, x);
matlabFunction(polar2Cartesian_jac, 'vars', {x}, 'file', 'polar2Cartesian_jac_fun.m');

%% reach error function
target = sym('target', [2, 1]);
reach_error = norm(polar2Cartesian - target)^2;
matlabFunction(reach_error, 'vars', {x, target}, 'file', 'reach_error_fun.m');
reach_error_jac = jacobian(reach_error, x);
matlabFunction(reach_error_jac, 'vars', {x, target}, 'file', 'reach_error_jac_fun.m');

%% 3-link signaling and reach error function
% this script generates: the amount of signaling function of 3-link arm
%                        the reach error in Cartesian
% it is seperated from the main dynamic generating code because the number 
% of states are different

% system dimension
n = 3;

% states and length
clear x l;
x = sym('x', [4*n, 1]);
l = sym('l', [3, 1]);

% position of joints
joints_pos = sym(zeros(4,3));
joints_pos(:, 1) = l(1)*[cos(x(1)); sin(x(1)); ...
                         cos(x(7)); sin(x(7))];
joints_pos(:, 2) = joints_pos(:, 1) + ...
                   l(2)*[cos(x(1)+x(2));sin(x(1)+x(2)); ...
                         cos(x(7)+x(8));sin(x(7)+x(8))];
joints_pos(:, 3) = joints_pos(:, 2) + ...
                   l(3)*[cos(x(1)+x(2)+x(3));sin(x(1)+x(2)+x(3)); ...
                         cos(x(7)+x(8)+x(9));sin(x(7)+x(8)+x(9))];
                     
% coordinate difference between corresponding joints
d = simplify([eye(2),-eye(2)]*joints_pos);

% weighted total distance of corresponding joints
signaling = simplify([1, 1, 1]*diag(transpose(d)*d));
matlabFunction(signaling, 'vars', {x,l}, 'file', 'signaling_fun.m');

% jacobian
signaling_jac = simplify(transpose(jacobian(signaling, x)));
matlabFunction(signaling_jac, 'vars', {x,l}, 'file', 'signaling_jac_fun.m');

% jacobian^2
signaling_jac2 = simplify(jacobian(signaling_jac, x));
matlabFunction(signaling_jac2, 'vars', {x,l}, 'file', 'signaling_jac2_fun.m');

% reach error 
target = sym('target', [4, 1]);
reach_error2 = norm(polar2Cartesian_fun(x(1:2*n)) - target(1:2))^2 + ...
               norm(polar2Cartesian_fun(x((2*n+1):4*n)) - target(3:4))^2;
matlabFunction(reach_error2, 'vars', {x, target}, 'file', 'reach_error2_fun.m');
reach_error2_jac = simplify(transpose(jacobian(reach_error2, x)));
matlabFunction(reach_error2_jac, 'vars', {x, target}, 'file', 'reach_error2_jac_fun.m');
reach_error2_jac2 = jacobian(reach_error2_jac, x);
matlabFunction(reach_error2_jac2, 'vars', {x, target}, 'file', 'reach_error2_jac2_fun.m');