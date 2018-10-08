%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project2_1                                             %%%
%%% Yiting Wang                                                      %%%
%%% Example 2.3 Page 31                                              %%%
%%% Target:                                                          %%%
%%%       Use the initialization and approching formula to calculate %%%
%%%       the value of partial equation.                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Initial number
dx = 0.1;
dt = 0.0025;
r = dt/dx^2;
x = 0:dx:1;
t = 0:dt:1;
M = length(t);
N = length(x) + 2; % Create 2 columns of the matrix
U = zeros(M,N);
U(1,:) = ones(1,N);

% Begin to generate the matrix
for i = 2:M
    U(i - 1,1) = U(i - 1,3) - 2 * dx * U(i - 1,2); % Calculate the beginning point
    U(i - 1,end) = U(i - 1,end - 2) - 2 * dx * U(i - 1,end - 1); % Calculate the ending point
    for j = 2:N-1
        U(i,j) = r * U(i - 1,j - 1) + (1 - 2 * r) * U(i - 1,j) + r * U(i - 1,j+1); % Follow the formula
    end
end
U(M,1) = U(M,3) - 2 * dx * U(M,2); % Calculate the beginning point when i = M
U(M,end) = U(M,end - 2) - 2 * dx * U(M,end - 1); % Calculate the ending point when i = M

%%% Output
fprintf('                                Table 2.10                   \n');
fprintf('------------------------------------------------------------------------------------\n');
fprintf('            i=0        i=1        i=2        i=3        i=4        i=5\n');
fprintf('            x=0        0.1        0.2        0.3        0.4        0.5\n');
fprintf('------------------------------------------------------------------------------------\n');
for k =0:M-1
    fprintf('t=%.4f | %.4f     %.4f     %.4f     %.4f     %.4f     %.4f\n', k*dt,U(k+1,2),U(k+1,3),U(k+1,4),U(k+1,5),U(k+1,6),U(k+1,7));
end

% Calculate the true value of the function
% Find about 1000 roots 
% Set the first points at 1, the plus pi*k(k = 1,2,3,...,1000)
% Use fzero to find the roots
% xx = [0:0.1:pi/2];
% yy1 = tan(xx);
% yy2 = 1/2./xx;
% plot(xx,yy1,'b.-',xx,yy2,'r--');
Max_points = 300;
begin_points = 1;
% Find points of alpha
for k = 1:Max_points
    alpha(k) = findPoints('x*tan(x)-1/2', 0 + (k - 1) * pi, pi/2 - 0.01 + (k - 1) * pi);
end
% Calculate the analytical result
u = ones(M,N-2);
for i = 2:M
    for j = 1:N-2
        sum = 0;
        for k=1:Max_points
            sum = sum + (sec(alpha(k)) / (3 + 4 * alpha(k)^2) * exp(- 4 * alpha(k)^2 * dt *(i - 1)) * cos(2 * alpha(k) * (dx * (j-1) -1/2)));
        end
        u(i,j) = 4 * sum;
    end
end
    
% Calculate the errors of numerical result and the analytical result
error = U(:,2:end-1) - u(:,:);
[X,T] = meshgrid(x,t);
mesh(x,t,error)
xlabel('x');
ylabel('t');
zlabel('error')
