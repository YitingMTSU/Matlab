%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 6                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the explict finite-difference formula                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% Initial number


xmin = 0;
dx = 0.1;
xmax = 1;

tmin = 0;
dt = 0.1;
tmax = 0.5;

x = xmin:dx:xmax; 
t = tmin:dt:tmax; 

M = length(t); % size of t
N = length(x); % size of x
u = zeros(N,M);

u(:,1) = 1/8*sin(pi*x);
h = dx;
k = dt;
r = k/h;

% Numerical solution
for j = 1:M-1
    if j == 1
        for i = 2:N-1
            u(i,j+1) = 1/2*(u(i-1,j)+u(i+1,j));
        end
    else
        for i = 2:N-1
            u(i,j+1) = r^2*u(i-1,j) + 2*(1-r^2)*u(i,j) + r^2*u(i+1,j) - u(i,j-1);
        end
    end
end

% print the output
fprintf('   x =      0       0.1       0.2       0.3       0.4       0.5\n');
fprintf('----------------------------------------------------------------\n');
for k = 2:M
    fprintf('t = %.1f | %.4f    %.4f    %.4f    %.4f    %.4f    %.4f\n',0.1*k,u(1,k),u(2,k),u(3,k),u(4,k),u(5,k),u(6,k));
end

% Analytical solution
U = zeros(N,M);
for j = 1:N
    for i = 1:M
        U(j,i) = 1/8*sin(pi*x(j))*cos(pi*t(i));
    end
end


[T,X] = meshgrid(t,x);
mesh(T,X,u);
xlabel('t');
ylabel('x');
zlabel('Numerical result');
title('Explict finite-difference formula');

figure

% plot the error at t=0.3
plot(x,u(:,4),'b-');
hold on
plot(x,U(:,4),'r-.');
xlabel('x');
ylabel('y');
title('Compare y-value between numerical and analytical solution');
legend('Numerical','Analytical');
figure
plot(x,u(:,4)-U(:,4),'g-');
xlabel('x');
ylabel('y');
title('The error at point t = 0.3');
        
