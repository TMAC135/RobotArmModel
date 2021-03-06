% This function returns linearized model of the system
% inputs: value of states: x
%         input to the system: u
% use the following functions generated by Symb_lagrange_solver.m
% M_fun.m
% C_fun.m
% phi_fun.m
% fric_fun.m
function [ A,B ] = linearized_model( x,u )
    
    % get the dimension of the system
    n = length(x)/2;
    % q = x(1:n);
    qdot = x((n+1):(2*n));
    
    % intermidiate steps
    M = M_fun(x);
    C = C_fun(x);
    phi = phi_fun([x;u]);
    fric = fric_fun(x);
    M_jac = M_jac_fun(x);
    C_jac = C_jac_fun(x);
    phi_jac_x = phi_jac_fun_x([x;u]);
    phi_jac_u = phi_jac_fun_u([x;u]);
    fric_jac = fric_jac_fun(x);
    
    % intermidiate steps
    invM_jac = zeros(n,n,2*n);
    for i=1:2*n
        invM_jac(:,:,i)=M\M_jac(:,:,i)/M;
    end
    
    U = -C*qdot - phi  - fric;
    
    qdot_jac = [zeros(n),eye(n)];
    
    %Calculate A and B
    A=[zeros(n),eye(n);zeros(n,2*n)];
    for i=1:2*n
        A(n+1:2*n,i) = -invM_jac(:,:,i)*U + ...
                       M\(-C_jac(:,:,i)*qdot - C*qdot_jac(:,i) - ...
                          phi_jac_x(:,i) - fric_jac(:,i));
    end
    
    B=zeros(2*n);
    for i=1:2*n
    B(n+1:2*n,i) = -M\phi_jac_u(:,i);
    end
end