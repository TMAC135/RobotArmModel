%run the linearlize model
function [X,U]=mpc_control(x_initial,x_target,u_initial,Q,R,N,dt)
%x_initial is the initial state of the function
%assuming the first half of the state is the location second half is the
%velocity
%x_target is the target of the function
% Q and R are cost matrices defined by x'*Q*x+u'*R*u
% You need to wirte in the form of [x(3) x(4) qddotfun(x)] for two link arm
% or [x(0)] for pendulum
% N is the time horizon 
% dt is the time interval
% l is the vector of the length
nx=size(x_initial,1);
nu=size(u_initial,1);
X=zeros(nx,N);%"recorder"
%initialize the initial state to the "recorder"
for i=1:nx
X(i,1)=x_initial(i);
end
U=zeros(nu,N);
for i=1:nu
U(i,1)=u_initial(i);
end
%The iteration of the LQR control
x=x_initial;
u=u_initial;
 for t=1:N-1
     [A1,B1]=linearized_model( x,u );
     A=eye(nx)+dt*A1;
     B=dt*B1;
     g=dt*(f(x,u)-A1*x-B1*u);
     u=LQR( nx,nu,A, B, Q, R,g,x,x_target);
     x=A*x+B*u+g;
     for i=1:nx
     X(i,t+1)=x(i);
     end
     for i=1:nu
     U(i,t+1)=u(i);
     end
     
 end
end