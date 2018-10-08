%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 4.1                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the explict, implict, Crank-Nicolson method            %%%
%%%       solve the patial equation and compare the error            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% fixed number in dx
dx = 0.01;
dt = 0.1;
xmax = 1;
tmax = 1;
max_step = 10;

U_EX_error = zeros(max_step,6);
U_IM_error = zeros(max_step,6);
U_CN_error = zeros(max_step,6);

for n = 1:max_step
    %Initial number
    dt = dt/ 2^(i-1);
    r = dt/dx^2;
    x = 0:dx:xmax;
    t = 0:dt:tmax;
    M = length(t);
    N = length(x) + 2;
    
    %-----------------------------------------------------------------------------------------------------
    % Explict method
    U_EX = zeros(M,N);
    U_EX(1,:) = ones(1,N);
     
    for i = 2:M
        U_EX(i-1,1) = U_EX(i-1,3) - 2*dx*U_EX(i-1,2); % Calculate the beginning point
        U_EX(i-1,end) = U_EX(i-1,end-2) - 2*dx*U_EX(i-1,end-1); % Calculate the ending point
        for j = 2:N-1
            U_EX(i,j) = r*U_EX(i-1,j-1) + (1-2*r)*U_EX(i-1,j) + r*U_EX(i-1,j+1); % Follow the formula
        end
    end
    U_EX(M,1) = U_EX(M,3) - 2*dx*U_EX(M,2); % Calculate the beginning point when i = M
    U_EX(M,end) = U_EX(M,end-2) - 2*dx*U_EX(M,end-1); % Calculate the ending point when i = M
    
    %-------------------------------------------------------------------------------------------
    % acturial answers
    Max_points = 300;
    % Find points of alpha
    for k = 1:Max_points
        %alpha(k) = findPoints('x*tan(x)-1/2', 0 + (k - 1) * pi, pi/2 - 0.01 + (k - 1) * pi);
        alpha(k) = fsolve('x*tan(x)-1/2',0.1+(k-1)*pi-eps);
    end
    % Calculate the analytical result
    u_ac = ones(M,N-2);
    for i = 2:M
        for j = 1:N-2
            sum = 0;
            for k=1:Max_points
                sum = sum + (sec(alpha(k)) / (3 + 4 * alpha(k)^2) * exp(- 4 * alpha(k)^2 * dt *(i - 1)) * cos(2 * alpha(k) * (dx * (j-1) -1/2)));
            end
            u_ac(i,j) = 4 * sum;
        end
    end

    %---------------------------------------------------------------------------------------------
    %Implict method
    
    U_IM = zeros(M,N);
    U_IM(1,:) = ones(1,N);

    % Construc matrix A
    A = zeros(N-2);
    for i=1:N-2
        if i==1
            A(i,i) = 1+2*r;
            A(i,i+1) = -r;
        elseif i==N-2
            A(i,i-1) = -r;
            A(i,i) = 1+2*r;
        else
            A(i,i-1) = -r;
            A(i,i) = 1+2*r;
            A(i,i+1) = -r;
        end
    end
    
    b = zeros(N-2,1);
    for i=2:M
        U_IM(i-1,1) = U_IM(i-1,3) - 2*dx*U_IM(i-1,2); % Calculate the beginning point
        U_IM(i-1,end) = U_IM(i-1,end-2) - 2*dx*U_IM(i-1,end-1); % Calculate the ending point
        b(1) = r*U_IM(i-1,1);
        b(end) = r*U_IM(i-1,end);
        u_temp = A\(transpose(U_IM(i-1,2:N-1))+b);
        U_IM(i,2:N-1) = transpose(u_temp);
    end
    U_IM(M,1) = U_IM(M,3) - 2*dx*U_IM(M,2); % Calculate the beginning point when i = M
    U_IM(M,end) = U_IM(M,end-2) - 2*dx*U_IM(M,end-1); % Calculate the ending point when i = M

    %----------------------------------------------------------------------------------------------------------%

    %%% Crank-Nicolson
    U_CN = zeros(M,N);
    U_CN(1,:) = ones(1,N);

    A1 = zeros(N - 2);
    A2 = zeros(N - 2);
    b1 = zeros(N-2,1);

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

    % Loop through time
    for j = 2:M
        U_CN(j-1,1) = U_CN(j-1,3) - 2*dx*U_CN(j-1,2); % Calculate the beginning point
        U_CN(j-1,end) = U_CN(j-1,end-2) - 2*dx*U_CN(j-1,end-1); % Calculate the ending point
        b1(1) = r*U_CN(j-1,1) + r*U_CN(j,1);
        b1(end) = r*U_CN(j-1,end) + r*U_CN(j,end);
        u_cntemp = A1\(A2*transpose(U_CN(j-1,2:N-1)) + b1);
        U_CN(j,2:N-1) = transpose(u_cntemp);
    end
    U_CN(M,1) = U_CN(M,3) - 2*dx*U_CN(M,2); % Calculate the beginning point when i = M
    U_CN(M,end) = U_CN(M,end-2) - 2*dx*U_CN(M,end-1); % Calculate the ending point when i = M
    
    %-----------------------------------------------------------------------------------------
    %calculate the error
    U_EX_ERROR = U_EX(:,2:N-1) - u_ac;
    U_IM_ERROR = U_IM(:,2:N-1) - u_ac;
    U_CN_ERROR = U_CN(:,2:N-1) - u_ac;
    
    U_EX_error(n,1) = sum(sum(abs(U_EX_ERROR)));
    U_EX_error(n,3) = sqrt(sum(sum(U_EX_ERROR.^2)));
    U_EX_error(n,5) = max(max(abs(U_EX_ERROR)));
    
    U_IM_error(n,1) = sum(sum(abs(U_IM_ERROR)));
    U_IM_error(n,3) = sqrt(sum(sum(U_IM_ERROR.^2)));
    U_IM_error(n,5) = max(max(abs(U_IM_ERROR)));
    
    U_CN_error(n,1) = sum(sum(abs(U_CN_ERROR)));
    U_CN_error(n,3) = sqrt(sum(sum(U_CN_ERROR.^2)));
    U_CN_error(n,5) = max(max(abs(U_CN_ERROR)));
    
end


for i=2:max_step
    U_EX_error(i,2) = log(U_EX_error(i,1)/U_EX_error(i-1,1))/log(1/2);
    U_IM_error(i,2) = log(U_IM_error(i,1)/U_IM_error(i-1,1))/log(1/2);
    U_CN_error(i,2) = log(U_CN_error(i,1)/U_CN_error(i-1,1))/log(1/2);
end
  