function u=LQR( nx,nu,A, B, Q, R,g,x,xd)
n=5;%the length of the mpc
A1=zeros(nx,nx,n);
B1=zeros(nx,nu,n);
Q1=zeros(nx,nx,n);
R1=zeros(nu,nu,n);
g1=zeros(nx,n);
l=zeros(nu,n);
h=zeros(nx,n);
for i=1:n
A1(:,:,i)=A;
B1(:,:,i)=B;
Q1(:,:,i)=Q;
R1(:,:,i)=R;
h(:,i)=-xd'*Q;
g1(:,i)=g;
end

[ K, s, ~, ~, ~ ] = LQR_design( A1, B1, Q1, R1, l, h, g1);


u=K(:,:,1)*x+s(:,1);
end