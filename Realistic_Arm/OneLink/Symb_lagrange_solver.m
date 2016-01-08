% this is a function to solve Euler-Lagrange Equations in a symbolic way
% 
% Two conditions are hold for this function: 
% 1. Kinetic energy is in a quadratic form of state qdot
% 2. Potential energy is independent of qdot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funtion input:
% Lagrangian function (L):a symbolic function of states q and qdot
%  
% State (q): a column vector, multi-dimensional, say the dimension is n
%  
% State (qdot): a column vector, represents the derivative of state q, has 
%               the same dimension as q
% 
% Torque(tor): a column vector, represents the input of the overall input, 
%              has the same dimension as q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function output:
% Mass matrix (M): a n-by-n matrix, same as D matrix in Tyler and Jinsun's 
%                  work
% 
% C matrix (C): a n-by-n matrix. C*qdot represents quadratic terms in qdot
%               in Eular-Lagrange equation
% 
% Terms concidering the change of potential energy (phi): a n-dimensional
%                                                         column vector
% 
% State dirivative functions (V): a 2n column vector shows the dirivative 
%                                 of q (which is qdot), and the dirivative 
%                                 of qdot (which is qddot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other notation: 
% Kinetic energy of the system (K): a 1-by-1 matrix
% Potential energy of the system(P): a 1-by-1 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code below is a 2-link arm dynamic testing example. aftering typing these
% lines, one should get matrix M, C, phi, and V.
% (Ctrl+R and Ctrl+T can help you comment or de-comment some lines of code)
%
% syms q qdot L tor;
% syms q1 q2 qdot1 qdot2 real;
% syms g l1 l2 m1 m2 tor1 tor2 I1 I2 real;
% q=[q1;q2];
% qdot=[qdot1;qdot2];
% tor=[tor1;tor2];
% L = 1/8*m1*l1^2*qdot1^2+...
%     .5*m2*(l1^2*qdot1^2+l1*l2*qdot1*(qdot1+qdot2)*cos(q2)+1/4*l2^2*(qdot1+qdot2)^2)+...
%     .5*I1*qdot1^2+.5*I2*(qdot1+qdot2)^2+...
%     m1*g*.5*l1*sin(q1)+m2*g*(l1*sin(q1)+.5*l2*sin(q1+q2));
% [M,C,phi,V]=Symb_lagrange_solver(L,q,qdot,tor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-link test
% clear all
% clc
% syms q qdot L tor;
% syms q1 q2 q3 qdot1 qdot2 qdot3 real;
% syms g l1 l2 m1 m2 tor1 tor2 tor3 I1 I2 real;
% q=[q1;q2;q3];
% qdot=[qdot1;qdot2;qdot3];
% tor=[tor1;tor2;tor3];
% L = 1*qdot1^2 + 5*qdot1*qdot2 + 3*qdot2^2 + 7*qdot3^2;
% [M,C,phi,V]=Symb_lagrange_solver(L,q,qdot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,C,phi] = Symb_lagrange_solver(T,V,x)

% get the system dimention
n = length(x)/2;
q = x(1:n);
qdot = x((n+1):(2*n));

%Calculate matrix M
M=simplify(jacobian(T,x((n+1):(2*n))));
M=simplify(jacobian(M,x((n+1):(2*n))));
n=max(size(M));%find the dimension of M

%Calculate matrix C

c = sym(zeros(n,n,n)); % initialization
C = sym(zeros(n,n));

for k=1:n
    for i=1:n
        for j=1:n
            c(i,j,k)=.5*(diff(M(k,j),q(i))+diff(M(k,i),q(j))-diff(M(i,j),q(k)));
        end
    end
end
for k=1:n
    for j=1:n
        for i=1:n
            C(k,j)=C(k,j)+c(i,j,k)*qdot(i);
        end
    end
end

C=simplify(C);

%Calculate jacobian of potential energy
phi=transpose(simplify(jacobian(V,q)));

% friction
% fric;

% save the expressions of M, C, phi and friction as functions
matlabFunction(M,    'vars',{x},'file','M_fun');
matlabFunction(C,    'vars',{x},'file','C_fun');
matlabFunction(phi,  'vars',{x},'file','phi_fun');
%%matlabFunction(fric, 'vars',{x},'file','fric_fun');

% functions for linearized model
% Jacobian of M w.r.t. x
M_jac = sym(zeros(n,n,2*n));
for i=1:2*n
     M_jac(:,:,i) = diff(M,x(i));
end
matlabFunction(M_jac,'vars',{x},'file','M_jac_fun');

% Jacobian of C w.r.t. x
C_jac = sym(zeros(n,n,2*n));
for i=1:2*n
    C_jac(:,:,i) = diff(C,x(i));
end
matlabFunction(C_jac,'vars',{x},'file','C_jac_fun');

% Jacobian of phi w.r.t. x
phi_jac = jacobian(phi,x);
matlabFunction(phi_jac,'vars',{x},'file','phi_jac_fun');

% Jacobian of friction w.r.t x
%%fric_jac = jacobian(fric,x);
%%matlabFunction(fric_jac,'vars',{x},'file','fric_jac_fun');

end