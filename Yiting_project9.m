%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 9                                              %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       Use the extrapolation method to solve PDF                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
clear;
clc;
 
%Initial set
xmin = 0;
xmax = 2;
ymin = 0;
ymax = 2;
tmin = 0;
tmax = 1;
h = 0.1; % step size of x and y
dT = 0.1; %step size of t

ERRback = zeros(1,6);
ERRPR = zeros(1,6);
ERRNS = zeros(1,6);
for m = 1:2
    for n = 1:3
        % set step of each vector
        dx = h/2^(n-1);
        dy = dx;
        dt = dT/10^(m-1);
        % generate the vector
        x = xmin:dx:xmax;
        y = ymin:dy:ymax;
        t = tmin:dt:tmax;
        % size of each vector
        nx = length(x)-2;
        ny = length(y)-2;
        nt = length(t)-1;
        
        ex = ones(nx,1);
        ey = ones(ny,1);
        Ax = diag(ex(1:end-1),-1)+diag(ex(1:end-1),1)+diag(-2*ex,0);
        Ay = diag(ey(1:end-1),-1)+diag(ey(1:end-1),1)+diag(-2*ey,0);
        
        Ix = eye(nx);
        Iy = eye(ny);
        
        A1 = kron(Ax,Iy)/dx^2;
        A2 = kron(Ix,Ay)/dy^2;
        
        %define the initial problem
        ind = symrcm(A1);
        u = ones(nx,1)*sin(pi*y(2:end-1)/2);
        uvec = u(:);
        
        
        % analytical result for t = tmax
        U = zeros(nx+2,ny+2);
        for i = 1:nx+2
            for j = 1:ny+2
                for k = 1:400 % take the sum of first 400 items
                    U(i,j) = U(i,j) + (1-(-1)^k)*4/k/pi*sin(k*pi*x(i)/2)*exp(-pi^2/4*(k^2+1)*tmax);
                end
                U(i,j) = sin(pi*y(j)/2)*U(i,j);
            end
        end
        
        % genetate the matrix P1 and P2
        P1 = (eye(size(A2)) - dt*A2);
        P2 = (eye(size(A1)) - dt*A1(ind,ind));
        [L1 U1] = lu(P1);
        [L2 U2] = lu(P2);
      
        
        % backward-difference method
        backvec = uvec; 
        for i = 1:nt
            backvec = L2\backvec(ind);
            backvec = U2\backvec;
            backvec = L1\backvec;
            backvec = U1\backvec;
        end
        backvec = backvec(ind(ind));
        Uback = U;
        Uback(2:end-1,2:end-1) = reshape(backvec,nx,ny);
        
        
        % Peaceman-Rachford method (FORMULA 3.7 PAGE 1219)
        % genetate the matrix P1 and P2
        P1 = (eye(size(A2)) - dt/2*A2);
        P2 = (eye(size(A1)) - dt/2*A1(ind,ind));
        [L1 U1] = lu(P1);
        [L2 U2] = lu(P2);
        PRvec = uvec;
        for i = 1:nt
            PRvec = L1\((eye(size(A1))+dt/2*A1)*PRvec);
            PRvec = U1\PRvec;
            PRvec = L2\((eye(size(A2))+dt/2*A2)*PRvec(ind));
            PRvec = U2\PRvec;
        end
        PRvec = PRvec(ind(ind));
        UPR = U;
        UPR(2:end-1,2:end-1) = reshape(PRvec,nx,ny);

        % Novel scheme (FORMULA 3.9-3.15 PAGE 1219)
        % genetate the matrix P1 P2 P3 and P4
        tao = dt/2;
        P1 = (eye(size(A2)) - tao*A2);
        P2 = (eye(size(A1)) - tao*A1(ind,ind));
        P3 = (eye(size(A2)) - 2*tao*A2);
        P4 = (eye(size(A1)) - 2*tao*A1(ind,ind));
        [L1 U1] = lu(P1);
        [L2 U2] = lu(P2);
        [L3 U3] = lu(P3);
        [L4 U4] = lu(P4);
        NSvec = uvec;
        for i = 1:nt
            NSvec0 = U1\(L1\(U2\(L2\NSvec(ind))));
            NSvec0 = U2\(L2\(U1\(L1\NSvec0))); % FORMULA 3.8
            NSvec1 = U4\(L4\(U3\(L3\NSvec(ind))));
            NSvec2 = U3\(L3\(U4\(L4\NSvec(ind))));
            NSvec = 2*NSvec0 - 1/2*(NSvec1 + NSvec2);
        end
        NSvec = NSvec(ind(ind));
        UNS = U;
        UNS(2:end-1,2:end-1) = reshape(NSvec,nx,ny);
        
        %find the error between approximate mehthod and true value
        ERRback(3*(m-1)+n) = max(max(abs(Uback-U)));
        ERRPR(3*(m-1)+n) = max(max(abs(UPR-U)));
        ERRNS(3*(m-1)+n) = max(max(abs(UNS-U)));
        % plot the result when dx=0.025, dt=0.1
        if dx == 0.025 && dt == 0.1
            [X Y] = meshgrid(x,y);
            mesh(X,Y,Uback);
            figure;
            mesh(X,Y,UPR);
            figure;
            mesh(X,Y,UNS);
        end
    end
end
fprintf('                                           Table 2\n');
fprintf('                     Errors in solving the two-dimension model problem at t=1.\n');
fprintf('---------------------------------------------------------------------------------------------------------\n');
fprintf('                        |              dt = 0.1              |              dt = 0.01             \n');
fprintf('         Method         | \th=0.1\t\th=0.05\t\th=0.025\t | \t h=0.1\t\t h=0.05\t\th=0.025\t     \n');
fprintf('---------------------------------------------------------------------------------------------------------\n');
fprintf('   Backward-difference  |  %.3e   %.3e   %.3e |  %.3e   %.3e   %.3e  \n',ERRback(1),...
    ERRback(2),ERRback(3),ERRback(4),ERRback(5),ERRback(6));
fprintf('   Peaceman-Rachford    |  %.3e   %.3e   %.3e |  %.3e   %.3e   %.3e  \n',ERRPR(1),...
    ERRPR(2),ERRPR(3),ERRPR(4),ERRPR(5),ERRPR(6));
fprintf('   Novel scheme         |  %.3e   %.3e   %.3e |  %.3e   %.3e   %.3e  \n',ERRNS(1),...
    ERRNS(2),ERRNS(3),ERRNS(4),ERRNS(5),ERRNS(6));
fprintf('---------------------------------------------------------------------------------------------------------\n');