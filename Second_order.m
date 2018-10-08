function [U] = Second_order(alpha,x,h,k,t,theta)

%%%%%%%%%%%%%%%%%%%%%%%%%  INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial set
beta = 0.8;

M = length(t) - 1; % size of t
N = length(x) - 1; % size of x
U = zeros(N + 1,M + 1);

gamma = 1/4 - theta/2;
I = eye(N-1);
%define the analytical function
syms X T
u = symfun(sqrt(2/(2+alpha*T))*exp(-(X-2-0.8*T)^2/(4*(2+alpha*T))),[X T]);
du = symfun(diff(u,T),[X T]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   CONSTRUCT MATRIX U  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N+1
    U(i,1) = u(x(i),0);
end
for j = 1:M+1
    U(1,j) = u(0,t(j));
    U(N+1,j) = u(x(end),t(j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   CONSTRUCT MATRIX du  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N+1
    dU(i,1) = du(x(i),0);
end
for j = 1:M+1
    dU(1,j) = du(0,t(j));
    dU(N+1,j) = du(x(end),t(j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   CONSTRUCT MATRIX A  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 1/h^2 * full(gallery('tridiag',N-1,alpha + beta*h/2, -2*alpha, alpha - beta*h/2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   CONSTRUCT VECTOR d  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta0 = 1/2 + theta;
beta1 = 1/2 - theta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   VECTOR b AND b'(db) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(N-1,M+1);
db = zeros(N-1,M+1);
b(1,:) = 1/h^2*(alpha + beta*h/2)*U(1,:);
b(end,:) = 1/h^2*(alpha - beta*h/2)*U(N+1,:);

db(1,:) = 1/h^2*(alpha + beta*h/2)*dU(1,:);
db(end,:) = 1/h^2*(alpha - beta*h/2)*dU(N+1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   SOLVE THE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:M
    % STEP 1 VECTOR d
    d = (I + k*beta0*A)*U(2:N,i) + k*beta0*b(:,i) + k*(beta1*I + k*theta*A)*b(:,i+1) + k^2*theta*db(:,i+1);
    Utemp = inv(I - gamma*k*A)*d;
    U(2:N,i+1) = inv(I - gamma*k*A)*Utemp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%