%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 3                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the explict, implict, Crank-Nicolson method            %%%
%%%       solve the patial equation.                                 %%%
%%% Equation:                                                        %%%
%%%          du/dt = d^2u/dx^2;                                      %%%
%%% Analvitical solution:                                            %%%
%%%          u|x=0 = 0 u|x=1 = 0                                     %%%
%%%          u|t=0 = 1                                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% Initial number
dx = 0.1;
dt = 0.001;
xmax = 1;
tmax = 1;
r = dt/dx^2;
x = 0:dx:xmax; % xmax 1
t = 0:dt:tmax; % tmax 1
M = length(t); % size of t
N = length(x); % size of x
U = zeros(M,N);
U(1,:) = ones(1,N);

%%% Explict method

u_ex = zeros(M,N);
u_ex(1,:) = ones(1,N);
for i = 2:M  % Calculate until j = 20
    for j = 2:N-1
        u_ex(i,j) = r*u_ex(i-1,j-1) + (1-2*r)*u_ex(i-1,j) + r*u_ex(i-1,j+1); % Follow the formula
    end
end

%%% Output
tableWith = 40;
for i=1:tableWith
    fprintf(' ');
end
fprintf('Table For Explict Method\n');
for i=1:tableWith*3
    fprintf('-');
end
fprintf('\n');
blankwith = 17; 
for i=1:blankwith
    fprintf(' ');
end
for i=1:N
    fprintf('i=%d        ',i-1);
end
fprintf('\n');
for i=1:blankwith
    fprintf(' ');
end
for i=1:N
    fprintf('x=%.1f      ',(i-1)*dx);
end
fprintf('\n');
for i=1:tableWith*3
    fprintf('-');
end
fprintf('\n');
for k =0:M-1
    fprintf('(j = %4d)t=%.3f|',k,k*dt);
    for i=1:N
        fprintf('%.4f     ',u_ex(k+1,i));
    end
    fprintf('\n');
end
[X,T] = meshgrid(x,t);
mesh(X,T,u_ex)
xlabel('x');
ylabel('t');
zlabel('Explict result')
title('Explict method')

%-------------------------------------------------------------------------------------------------------%

%%% Implict method

u_im = zeros(M,N);
u_im(1,:) = ones(1,N);

% Construc matrix A
A = zeros(N);
for i=1:N
    if i==1
        A(i,i) = 1+2*r;
        A(i,i+1) = -r;
    elseif i==N
        A(i,i-1) = -r;
        A(i,i) = 1+2*r;
    else
        A(i,i-1) = -r;
        A(i,i) = 1+2*r;
        A(i,i+1) = -r;
    end
end

for i=2:M
    u_temp = A\transpose(u_im(i-1,:));
    u_im(i,:) = transpose(u_temp);
    u_im(i,1) = 0;
    u_im(i,end) = 0;
end

figure(2)
mesh(X,T,u_im)
xlabel('x');
ylabel('t');
zlabel('Implict result')
title('Implict method')



%----------------------------------------------------------------------------------------------------------%

%%% Crank-Nicolson

u_cn = zeros(M,N);
u_cn(1,:) = ones(1,N);

A1 = zeros(N - 2);
A2 = zeros(N - 2);
b = zeros(N-2,1);

% A1 defination
A1(1,1) = 2 + 2*r;
A1(1,2) = -r;
for i = 2:N-3
    A1(i,i-1) = -r;
    A1(i,i) = 2 + 2*r;
    A1(i,i+1) = -r;
end
A1(N-2,N-3) = -r;
A1(N-2,N-2) = 2 + 2*r;

% A2 defination
A2(1,1) = 2 - 2*r;
A2(1,2) = r;
for i = 2:N-3
    A2(i,i-1) = r;
    A2(i,i) = 2 - 2*r;
    A2(i,i+1) = r;
end
A2(N-2,N-3) = r;
A2(N-2,N-2) = 2 - 2*r;

% b defination


% Loop through time
for j = 2:M
    b(1) = r*u_cn(j-1,1) + r*u_cn(j,1);
    b(end) = r*u_cn(j-1,end) + r*u_cn(j,end);
    u_cntemp = A1\(A2*transpose(u_cn(j-1,2:N-1)) + b);
    u_cn(j,2:N-1) = transpose(u_cntemp);
end

figure(3)
mesh(X,T,u_cn)
xlabel('x');
ylabel('t');
zlabel('Crank-Nicolson result')
title('Crank-Nicolson')