%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 7.2                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the explict finite-difference formula                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% Initial number
k = 0.1; 
ALPHA = [0.005 0.0005];
beta = 0.8;
H = [0.1,0.05,0.025,0.0125];
M = 10;

xmin = 0;
xmax = 1;
tmin = 0;

fprintf('                    Table II. Maximum absolute error at t = 1 using k = 0.1    \n');
fprintf('------------------------------------------------------------------------------------------------------------------\n');
fprintf('alpha                             0.005                                          0.0005\n');
fprintf('         -------------------------------------------------- ------------------------------------------------------\n');
fprintf('  h        BDF2         MDER         SLMM         SLMP        BDF2         MDER         SLMM         SLMP\n');
fprintf('------------------------------------------------------------------------------------------------------------------\n');
for i = 1:length(H)
    for j = 1:length(ALPHA)
        alpha = ALPHA(j);
        h = H(i);      
        dx = h;
        dt = k;
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
        
        % Matrix A1(BDF2)
        A1 = 1/h^2 * full(gallery('tridiag',N-1,alpha + beta*h/2, -2*alpha, alpha - beta*h/2));

        %vector b
        b = zeros(N-1,M);
        b(1,:) = 1/h^2*(alpha + beta*h/2)*u(1,2:end);
        b(end,:) = 1/h^2*(alpha - beta*h/2)*u(N+1,2:end);
        
        U_BDF2 = u;
        for l = 1:M
            U_BDF2(2:N,l+1) = inv(A1)*(du(:,l) - b(:,l));
        end
        
        err_BDF2 = max(abs(U_exact(:,end) - U_BDF2(:,end)));
        if h == 0.0125 && alpha == 0.005
            plot(x,U_exact(:,3),'ro');
            hold on
            plot(x,U_BDF2(:,3),'g-');
            title('Compare solution between BDF2 and analytical at t=0.3 when h = 0.0125, alpha = 0.005');
            legend('Analytical','BDF2');
            figure
        end
        %--------------------------------------------------------
        %MDER
        
        theta = -3/2-sqrt(2);
        U_MDER = Second_order(alpha,x,h,k,t,theta);
        err_MDER = max(abs(U_exact(:,end) - U_MDER(:,end)));
        if h == 0.0125 && alpha == 0.005
            plot(x,U_exact(:,3),'ro');
            hold on
            plot(x,U_MDER(:,3),'g-');
            title('Compare solution between MDER and analytical at t=0.3 when h = 0.0125, alpha = 0.005');
            legend('Analytical','MDER');
            figure
        end
        
        %---------------------------------------------------------
        %SLMM
        
        theta = (-1-sqrt(2))/2;
        U_SLMM = LMM(alpha,x,h,k,t,theta);
        err_SLMM = max(abs(U_exact(:,end) - U_SLMM(:,end)));
        if h == 0.0125 && alpha == 0.005
            plot(x,U_exact(:,3),'ro');
            hold on
            plot(x,U_SLMM(:,3),'g-');
            title('Compare solution between SLMM and analytical at t=0.3 when h = 0.0125, alpha = 0.005');
            legend('Analytical','SLMM');
            figure
        end
        
        %---------------------------------------------------------
        %SLMP
        
        theta = (-1+sqrt(2))/2;
        U_SLMP =LMM(alpha,x,h,k,t,theta);
        err_SLMP = max(abs(U_exact(:,end) - U_SLMP(:,end)));
        if h == 0.0125 && alpha == 0.005
            plot(x,U_exact(:,3),'ro');
            hold on
            plot(x,U_SLMP(:,3),'g-');
            title('Compare solution between SLMP and analytical at t=0.3 when h = 0.0125, alpha = 0.005');
            legend('Analytical','SLMP');
            figure
        end
        
        
        %print the result
        if j == 1
            fprintf(' %.4f   %.7f    %.7f    %.7f    %.7f',h,err_BDF2,err_MDER,err_SLMM,err_SLMP);
        else
            fprintf('   %.7f    %.7f    %.7f    %.7f\n',err_BDF2,err_MDER,err_SLMM,err_SLMP);
%             Rate_BDF2 = log(err_BDF2(2)/err_BDF2(1))/log(1/2);
%             Rate_BDF2w = log(err_BDF2w(2)/err_BDF2w(1))/log(1/2);
%             Rate_MDER = log(err_MDER(2)/err_MDER(1))/log(1/2);
%             Rate_SLMM = log(err_SLMM(2)/err_SLMM(1))/log(1/2);
%             Rate_SLMP = log(err_SLMP(2)/err_SLMP(1))/log(1/2);
%             fprintf(' %d   %.2f  %.3f   %.7f  %.4f  %.7f  %.4f  %.7f  %.4f  %.7f  %.4f  %.7f  %.4f \n',M,h,k,err_BDF2w(j),Rate_BDF2w,err_BDF2(j),Rate_BDF2,err_MDER(j),Rate_MDER,err_SLMM(j),Rate_SLMM,err_SLMP(j),Rate_SLMP);
        end
    end
    
    
end
fprintf('------------------------------------------------------------------------------------------------------------------\n')