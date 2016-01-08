%%The following code is to realize the reaching with 2 springs in
%%the single link. Without any input to the system, we regard the [k1;k2]
%%as our input(u) which are two spring coefficients. By giving control
%%laws to this virtual input, we are possible to control the arm. In this
%%code, we use the  generarized iterative LQR method and it does show us
%%great performance for this model. 

clear all;

%%the simple arm with two spring attached in the link
 syms q qdot;
 syms g l m I r r1 r2 k1 k2 d c s; %k1 and k2 are spring coefficient 
                                   %for the two springs
 

%% arm parameters
m = 1;
l = 1; 
c = l/2;  
I = 0.5.*m.*l.^2;
g=0;
d=1;
s=0.8*c;
%x_init = [pi/2; 0];
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
 fric=0;
%fric=0.5*qdot;
%P=matlabFunction(fric);
x=[q;qdot;k1;k2];%state variable for the system,
[M,C,phi] = Symb_lagrange_solver_varying_k(T,V,x);
phi_func=matlabFunction(phi);


%% *Simulation for targrt reaching*
dt = 1e-2; % sample time
N = 600; % total steps
t = 0:dt:((N-1)*dt);


% cost function parameters
nx=2;
nu=2;
x_init = [pi/2; 0];
x_target = [pi/1.2; 0];
u_init = [0; 0];

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

Q(:, :, N) = 10000*eye(nx, nx); % final cost for the target reaching
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
            A(:, :, k)=[0 ,1;0,0];
            phi=phi_func(k1,k2,pi/2);
            %phi_jac=jacobian(phi,[k1,k2]);
            B(:, :, k)=[0,0;-M\(jacobian(phi,[k1;k2]))];
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
            A(:, :, k)=[0 ,1;0,0];
            phi=phi_func(k1,k2,x(1,k,kk-1));
            B(:, :, k)=[0,0;-M\(jacobian(phi,[k1;k2]))];
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
    for kk=1:page
      D=threelink_vert2D_createmovie_reaching(x(:,:,kk),t,10,l,x_target);      
    end
    
        


