

function [ K, s, P, p, c ] = LQR_design( A, B, Q, R, l, h, g)
% LQR_DESIGN Summary of this function goes here

% cost funtion
%  N
% sum (ui^T*Ri*ui + 2*li*ui + xi^T*Qi*xi + 2*hi*xi) + xN^T*QN*xN + hN*xN
% i=1

% system dynamics
% x(i+1) = Ai*xi + Bi*ui

% Check parameters dimensions
nx = size(A, 1);
nu = size(B, 2);
N = size(A, 3);

% Initializations
P = zeros(nx, nx, N);
P(:, :, N) = Q(:, :, N);
p = zeros(nx, N);
p(:, N) = h(:, N);
K = zeros(nu, nx, N);
s = zeros(nu, N);
c = zeros(1, N);

AT = zeros(nx, nx, N);
BT = zeros(nu, nx, N);
for k = 1:N
    AT(:, :, k) = transpose(A(:, :, k));
    BT(:, :, k) = transpose(B(:, :, k));
end

% back propagation
for k = N:(-1):2
    if min(eig((R(:, :, k-1)+BT(:, :, k-1)*P(:, :, k)*B(:, :, k-1)))) < 1e-10
        c(1) = -inf;
        break;
    end
    P(:, :, k-1) = Q(:, :, k-1) + ...
        AT(:, :, k-1)*P(:, :, k)*A(:, :, k-1) - ...
        AT(:, :, k-1)*P(:, :, k)*B(:, :, k-1)/ ...
        (R(:, :, k-1)+BT(:, :, k-1)*P(:, :, k)*B(:, :, k-1))* ...
        BT(:, :, k-1)*P(:, :, k)*A(:, :, k-1);
    p(:, k-1) = h(:, k-1) + AT(:, :, k-1)*p(:, k) + ...
        AT(:, :, k-1)*P(:, :, k)*g(:, k-1)- ...
        AT(:, :, k-1)*P(:, :, k)*B(:, :, k-1)/ ...
        (R(:, :, k-1)+BT(:, :, k-1)*P(:, :, k)*B(:, :, k-1))* ...
        (l(:, k-1) + BT(:, :, k-1)*p(:, k) + ...
        BT(:,:,k-1) * P(:,:,k)* g(:,k-1));
    c(k-1) = c(k) - (transpose(l(:, k-1)) + transpose(p(:, k))*B(:, :, k-1) + ...
        transpose(g(:,k-1)) * P(:,:,k) * B(:,:,k-1))/ ...
        (R(:, :, k-1) + BT(:, :, k-1)*P(:, :, k)*B(:, :, k-1))* ...
        (l(:, k-1) + BT(:, :, k-1)*p(:, k) + ...
        BT(:,:,k-1) * P(:,:,k) * g(:,k-1)) + ...
        2*transpose(p(:, k))*g(:, k-1) + transpose(g(:, k-1))*P(:, :, k)*g(:, k-1);
end
% optimal control law
if c(1)>-inf
    for k = 1:(N-1)
        K(:, :, k) = - (R(:, :, k) + BT(:, :, k)*P(:, :, k+1)*B(:, :, k))\ ...
            BT(:, :, k)*P(:, :, k+1)*A(:, :, k);
        s(:, k) = - (R(:, :, k) + BT(:, :, k)*P(:, :, k+1)*B(:, :, k))\ ...
            (l(:, k) + BT(:, :, k)*p(:, k+1) + ...
            BT(:, :, k)*P(:, :, k+1)*g(:, k));
    end
end
end