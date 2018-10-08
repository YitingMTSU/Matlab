%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project1                                               %%%
%%% Example 2.1 Page 13                                              %%%
%%% Target:                                                          %%%
%%%       Use the initialization and approching formula to calculate %%%
%%%       the value of partial equation.                             %%%
%%% Equation:                                                        %%%
%%%          u(i,j+1) = 1/10 * (u(i-1,j) + 8 * u(i,j) + u(i+1,j));   %%%
%%% Analvitical solution:                                            %%%
%%%          u = 8/pi^2 *                                            %%%
%%%             sum(n)(1/n^2)(sin1/2*n*pi)(sinn*pi*x)exp(-n^2pi^2t); %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

x = [0:0.1:0.6];  % Define the range of x
Num = length(x);  % Calculate the size of x
u(1,:) = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8];  % Initialize the u when t = 0
for j = 1:20  % Calculate until j = 20
    for i = 1:Num
        if i == 1
            u(j+1,i) = 0; % u(x,t) = 0 when x = 0
        elseif i == Num
            u(j+1,i) = u(j+1,Num - 2); % u(x,t) is symmetric by x = 0.5
        else
            u(j+1,i) = 1/10 * (u(j,i-1) + 8 * u(j,i) + u(j,i+1)); % Follow the formula
        end
    end
end

%%% True value of the function when x = 0.3
xx = 0.3;
for j = 1:20  % Calculate until j = 20
    sum = 0;
    for n = 1:1000
        sum = sum + 1/n^2*sin(1/2*n*pi)*sin(n*pi*xx)*exp(-n^2*pi^2*j/1000);
    end
    U03(j) = sum*8 / pi^2;
end

%%% True value of the function when x = 0.5
xx = 0.5;
for j = 1:20  % Calculate until j = 20
    sum = 0;
    for n = 1:1000
        sum = sum + 1/n^2*sin(1/2*n*pi)*sin(n*pi*xx)*exp(-n^2*pi^2*j/1000);
    end
    U05(j) = sum*8 / pi^2;
end

%%% Output
fprintf('                                         Table 2.2                   \n');
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('                 i=0       i=1       i=2       i=3       i=4       i=5       i=6\n');
fprintf('                 x=0       0.1       0.2       0.3       0.4       0.5       0.6\n');
fprintf('--------------------------------------------------------------------------------------------\n');
for k =0:20
    fprintf('(j = %2d)t=%.3f |%.4f     %.4f     %.4f     %.4f     %.4f     %.4f     %.4f\n', k,k*0.001,u(k+1,1),u(k+1,2),u(k+1,3),u(k+1,4),u(k+1,5),u(k+1,6),u(k+1,7));
end