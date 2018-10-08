%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 4.2                                            %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the extrapolation method                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% Initial number
dx = 0.025;
dt = 0.025;
xmax = 1;
tmax = 1;
r1 = dt/dx^2;
r2 = r1/2;
x = 0:dx:xmax; % xmax 1
t = 0:dt:tmax; % tmax 1
M = length(t); % size of t
N = length(x); % size of x
u2 = zeros(M,N);
u2(1,2:N-1) = ones(1,N-2);
u1 = zeros(M,N);
u1(1,2:N-1) = ones(1,N-2);

A1 = zeros(N - 2);
A2 = zeros(N - 2);

I = eye(N-2);

% A1 defination
A1(1,1) = 1 + 4*r1;
A1(1,2) = -2*r1;
for i = 2:N-3
    A1(i,i-1) = -2*r1;
    A1(i,i) = 1 + 4*r1;
    A1(i,i+1) = -2*r1;
end
A1(N-2,N-3) = -2*r1;
A1(N-2,N-2) = 1 + 4*r1;

% A2 defination
A2(1,1) = 1 + 2*r2;
A2(1,2) = -r2;
for i = 2:N-3
    A2(i,i-1) = -r2;
    A2(i,i) = 1 + 2*r2;
    A2(i,i+1) = -r2;
end
A2(N-2,N-3) = -r2;
A2(N-2,N-2) = 1 + 2*r2;


% Loop through time
for j = 2:M
    u1_temp = A1\(transpose(u1(j-1,2:N-1)));
    u1(j,2:N-1) = transpose(u1_temp);
end

% Loop through time
for j = 2:M
    u2_temp = A2\(transpose(u2(j-1,2:N-1)));
    u2_temp = A2\u2_temp;
    u2(j,2:N-1) = transpose(u2_temp);
end


u_ac = 2*u2-u1;
[X,T] = meshgrid(x,t);
mesh(x,t,u_ac)
xlabel('x');
ylabel('t');
zlabel('extrapolation')
