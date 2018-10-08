%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project2_2                                             %%%
%%% Yiting Wang                                                      %%%
%%% Target:                                                          %%%
%%%       Use the matrix method to calculate                         %%%
%%%       the value of partial equation.                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Initialization
T = 1;
a = 1;
theta = 0;
dt = 0.1;
dx = 0.1;
k = 1;
r = k * dt / dx / dx;

% Spatial points
x = 0:dx:a;
% Temporal points
t = 0:dt:T;

N = length(x);
M = length(t);

U = zeros(N, M);

U(:,1) = sin(pi*x'/a);

A1 = zeros(N - 2);
A2 = zeros(N - 2);

% A1 defination
A1(1,1) = 1 + 2 * r * theta;
A1(1,2) = - r * theta;
for i = 2:N-3
    A1(i,i-1) = -r*theta;
    A1(i,i) = 1 + 2*r*theta;
    A1(i,i+1) = -r*theta;
end
A1(N-2,N-3) = -r*theta;
A1(N-2,N-2) = 1 + 2*r*theta;

% A2 defination
A2(1,1) = 1 - 2*r*(1-theta);
A2(1,2) = r * (1-theta);
for i = 2:N-3
    A2(i,i-1) = r*(1-theta);
    A2(i,i) = 1 - 2*r*(1-theta);
    A2(i,i+1) = r*(1-theta);
end
A2(N-2,N-3) = r*(1-theta);
A2(N-2,N-2) = 1 - 2*r*(1-theta);

% Loop through time
for j = 2:M
    U(2:N-1,j) = A1\(A2*U(2:N-1,j-1));
end

plot(x',U(:,1));
hold on
plot(x',U(:,end),'--o');

figure(2)
[T,X] = meshgrid(x,t);
mesh(T,X,U');