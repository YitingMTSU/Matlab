%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 7                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the explict finite-difference formula                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% Initial number
k = 0.05; 
alpha = 0.005;
beta = 0.8;

MM = [15 20];

xmin = 0;
xmax = 1;
tmin = 0;

fprintf('                    Table 1. Maximum absolute error and rate of convergence for alpha = 0.005   \n');
fprintf('------------------------------------------------------------------------------------------------------------------\n');
fprintf('alpha = 0.005            BDF2w              BDF2               MDER                SLMM               SLMP\n');
fprintf('-----------------  ------------------ ------------------ ------------------- ------------------ ------------------\n');
fprintf('  M     h     k     Error      Rate    Error      Rate    Error       Rate    Error      Rate    Error      Rate   \n');
fprintf('------------------------------------------------------------------------------------------------------------------\n');
for i = 1:length(MM)
    for j = 1:2
        k = k/j;
        h = 2*k;      
        dx = h;
        dt = k;
        M = MM(i)*j;
        x = xmin:dx:xmax;
        t = tmin:dt:dt*M; 
        N = length(x) - 1;
        u = zeros(N + 1,M + 1);
        %define the analytical function
        syms X T
        fu = symfun(sqrt(2/(2+alpha*T))*exp(-(X-2-0.8*T)^2/(4*(2+alpha*T))),[X T]);
        Dfu_DT = symfun(diff(fu,T),[X T]);
        % Initial  and boundary conditions
        for l = 1:N+1
            u(l,1) = fu(x(l),0);
        end
        for l = 1:M+1
            u(1,l) = fu(0,t(l));
            u(N+1,l) = fu(xmax,t(l));
        end
        %Derivative matrix
        for p = 1:N-1
            for q = 1:M
                du(p,q) = eval(Dfu_DT(x(p+1),t(q+1)));
            end
        end
        
        %-------------------------------------------------------------------------------
        %exact value
        U_exact = u;
        for p = 1:N+1
            for q = 1:M+1
                U_exact(p,q) = eval(fu(x(p),t(q)));
            end
        end
        %-------------------------------------------------------------------------------%
        %BDF2
        
        % Matrix A1(BDF2) A2(BDF2w)
        A1 = 1/h^2 * full(gallery('tridiag',N-1,alpha + beta*h/2, -2*alpha, alpha - beta*h/2));
        A2 = 1/h^2 * full(gallery('tridiag',N-1,alpha + beta*h,-(2*alpha + beta*h),alpha));

        %vector b
        b = zeros(N-1,M);
        b(1,:) = 1/h^2*(alpha + beta*h/2)*u(1,2:end);
        b(end,:) = 1/h^2*(alpha - beta*h/2)*u(N+1,2:end);
        
        U_BDF2 = u;
        for l = 1:M
            U_BDF2(2:N,l+1) = inv(A1)*(du(:,l) - b(:,l));
        end
        
        err_BDF2(j) = AbsErr(U_exact,U_BDF2);
        %--------------------------------------------------------
        %BDF2w
        %vector b
        b = zeros(N-1,M);
        b(1,:) = 1/h^2*(alpha + beta*h)*u(1,2:end);
        b(end,:) = 1/h^2*(alpha)*u(N+1,2:end);
        
        U_BDF2w = u;
        for l = 1:M
            U_BDF2w(2:N,l+1) = inv(A2)*(du(:,l) - b(:,l));
        end
        err_BDF2w(j) = AbsErr(U_exact,U_BDF2w);
        
        %--------------------------------------------------------
        %MDER
        
        theta = -3/2-sqrt(2);
        U_MDER = Second_order(alpha,x,h,k,t,theta);
        err_MDER(j) = AbsErr(U_exact,U_MDER);
        
        %---------------------------------------------------------
        %SLMM
        
        theta = (-1-sqrt(2))/2;
        U_SLMM = LMM(alpha,x,h,k,t,theta);
        err_SLMM(j) = AbsErr(U_exact,U_SLMM);
        
        %---------------------------------------------------------
        %SLMP
        
        theta = (-1+sqrt(2))/2;
        U_SLMP =LMM(alpha,x,h,k,t,theta);
        err_SLMP(j) = AbsErr(U_exact,U_SLMP);
        
        
        %print the result
        if j == 1
            fprintf(' %d   %.2f  %.3f   %.7f          %.7f          %.7f          %.7f          %.7f\n',M,h,k,err_BDF2w(j),err_BDF2(j),err_MDER(j),err_SLMM(j),err_SLMP(j));
        else
            Rate_BDF2 = log(err_BDF2(2)/err_BDF2(1))/log(1/2);
            Rate_BDF2w = log(err_BDF2w(2)/err_BDF2w(1))/log(1/2);
            Rate_MDER = log(err_MDER(2)/err_MDER(1))/log(1/2);
            Rate_SLMM = log(err_SLMM(2)/err_SLMM(1))/log(1/2);
            Rate_SLMP = log(err_SLMP(2)/err_SLMP(1))/log(1/2);
            fprintf(' %d   %.2f  %.3f   %.7f  %.4f  %.7f  %.4f  %.7f  %.4f  %.7f  %.4f  %.7f  %.4f \n',M,h,k,err_BDF2w(j),Rate_BDF2w,err_BDF2(j),Rate_BDF2,err_MDER(j),Rate_MDER,err_SLMM(j),Rate_SLMM,err_SLMP(j),Rate_SLMP);
        end
    end
    
    
end
fprintf('------------------------------------------------------------------------------------------------------------------\n')