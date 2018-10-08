%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMS 7300 Project 10                                             %%%
%%%                                                                  %%%
%%% Target:                                                          %%%
%%%       GLOBAL EXTRAPOLATION OF A FIRST ORDER SPLITING METHOD      %%% 
%%%       Use the extrapolation method to solve PDF FOR 3-D          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
clear;
clc;
 
%Initial set
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
zmin = 0;
zmax = 1;
tmin = 0;
tmax = 1;
gamma = 1/6;
h = 0.1; % step size of x and y
dT = 1/12; %step size of t

ERRback = zeros(1,2);
ERRLM = zeros(1,2);
ERRGE = zeros(1,2);
for m = 1:1
    for n = 1:1
        % set step of each vector
        dx = h/2^(n-1);
        dy = dx;
        dz = dx;
        dt = dT/10^(m-1);
        % generate the vector
        x = xmin:dx:xmax;
        y = ymin:dy:ymax;
        z = zmin:dy:zmax;
        t = tmin:dt:tmax;
        % size of each vector
        nx = length(x)-2;
        ny = length(y)-2;
        nz = length(z)-2;
        nt = length(t)-1;
        
        ex = ones(nx,1);
        ey = ones(ny,1);
        ez = ones(nz,1);
        Ax = diag(ex(1:end-1),-1)+diag(ex(1:end-1),1)+diag(-2*ex,0);
        Ay = diag(ey(1:end-1),-1)+diag(ey(1:end-1),1)+diag(-2*ey,0);
        Az = diag(ez(1:end-1),-1)+diag(ez(1:end-1),1)+diag(-2*ez,0);
        
        Ix = eye(nx);
        Iy = eye(ny);
        Iz = eye(nz);
        
        A1 = kron(kron(Ax,Iy),Iz)/dx^2*gamma;
        A2 = kron(kron(Ix,Ay),Iz)/dy^2*gamma;
        A3 = kron(kron(Ix,Iy),Az)/dz^2*gamma;
        
        %define the initial problem
        ind1 = symrcm(A1);
        ind2 = symrcm(A2);
        for i = 1:nx
            for j = 1:ny
                for k = 1:nz
                    u(i,j,k) = sin(pi*y(j+1))*sin(pi*z(k+1));
                end
            end
        end

        uvec = u(:);
        
        
        % analytical result for t = tmax
        U = zeros(nx+2,ny+2,nz+2);
        for i = 1:nx+2
            for j = 1:ny+2
                for k = 1:nz+2
                    for num = 1:400 % take the sum of first 400 items
                        U(i,j,k) = U(i,j,k) + (1-(-1)^num)*2/num/pi*sin(num*pi*x(i))*exp(-gamma*pi^2*(num^2+2)*tmax);
                    end
                    U(i,j,k) = sin(pi*y(j))*sin(pi*z(k))*U(i,j,k);
                end
            end
        end
        %plot the analytical solution
        [X Y Z] = meshgrid(x,y,z);
        
        for i = 1:nz+2
            mesh(X(:,:,i),Y(:,:,i),U(:,:,i))
            zlim([0, 0.013]);
            title('Analytical Solution at t = 1');
            pause(0.25);
        end
        figure;

      
        
        % backward difference method
        % genetate the matrix P1 and P2
        begin = cputime;
        P1 = eye(size(A3)) - dt*A3;
        P2 = eye(size(A2)) - dt*A2(ind2,ind2);
        P3 = eye(size(A1)) - dt*A1(ind1,ind1);

        backvec = uvec; 
        for i = 1:nt
            backvec = P3\backvec;
            backvec = P2\backvec;
            backvec = P1\backvec;
        end
        
        Uback = U;
        Uback(2:end-1,2:end-1,2:end-1) = reshape(backvec,nx,ny,nz);
        ERRback(2) = cputime - begin;
        % plot the result for backward-difference scheme 
        for i = 1:nz+2
            mesh(X(:,:,i),Y(:,:,i),Uback(:,:,i))
            zlim([0, 0.013]);
            title('Backward-Difference scheme at t = 1');
            pause(0.25);
        end
        figure;
        
        % Lawson-Morris scheme (FORMULA 4.3 PAGE 774)
        % genetate the matrix P1 and P2
        begin = cputime;
        tao = dt/2;
        P1 = eye(size(A3)) - tao*A3;
        P2 = eye(size(A2)) - tao*A2(ind2,ind2);
        P3 = eye(size(A1)) - tao*A1(ind1,ind1);
        P4 = eye(size(A3)) - 2*tao*A3;
        P5 = eye(size(A2)) - 2*tao*A2(ind2,ind2);
        P6 = eye(size(A1)) - 2*tao*A1(ind1,ind1);

        
        LMvec = uvec;
        for i = 1:nt
            LMvec0 = P1\(P2\(P3\(P3\(P2\(P1\LMvec)))));
            LMvec1 = P4\(P5\(P6\LMvec));
            LMvec2 = P6\(P5\(P4\LMvec));
            LMvec = 2*LMvec0 - 1/2*(LMvec1+LMvec2);
        end
        
        ULM = U;
        ULM(2:end-1,2:end-1,2:end-1) = reshape(LMvec,nx,ny,nz);
        ERRLM(2) = cputime - begin;
        
        % plot the result for Lawson-Morris scheme 
        for i = 1:nz+2
            mesh(X(:,:,i),Y(:,:,i),ULM(:,:,i))
            zlim([0, 0.013]);
            title('Lawson-Morris scheme at t = 1');
            pause(0.25);
        end
        figure;

        % Global extrap (FORMULA 3.9-3.15 PAGE 1219)
        
        begin = cputime;
        P1 = eye(size(A3)) - dt/3*A3;
        P2 = eye(size(A2)) - dt/3*A2(ind2,ind2);
        P3 = eye(size(A1)) - dt/3*A1(ind1,ind1);
        P4 = eye(size(A3)) - dt/2*A3;
        P5 = eye(size(A2)) - dt/2*A2(ind2,ind2);
        P6 = eye(size(A1)) - dt/2*A1(ind1,ind1);
        P7 = eye(size(A3)) - dt*A3;
        P8 = eye(size(A2)) - dt*A2(ind2,ind2);
        P9 = eye(size(A1)) - dt*A1(ind1,ind1);

        GEvec = uvec;
        for i = 1:nt
            GEvec1 = P7\(P8\(P9\GEvec));
            
            GEvec2 = P6\(P5\(P4\(P4\(P5\(P6\GEvec)))));
            
            GEvec3 = P1\(P2\(P3\(P3\(P2\(P1\(P1\(P2\(P3\GEvec))))))));
            
            GEvec = 0.5*GEvec1 - 4*GEvec2 + 4.5*GEvec3;
        end
        
        UGE = U;
        UGE(2:end-1,2:end-1,2:end-1) = reshape(GEvec,nx,ny,nz);
        ERRGE(2) = cputime - begin;
        
        % plot the result for Global extrap scheme 
        for i = 1:nz+2
            mesh(X(:,:,i),Y(:,:,i),ULM(:,:,i))
            zlim([0, 0.013]);
            title('Global extrap scheme at t = 1');
            pause(0.25);
        end
        
        %find the error between approximate mehthod and true value
        ERRback(3*(m-1)+n) = max(max(max(abs(Uback-U))));
        ERRLM(3*(m-1)+n) = max(max(max(abs(ULM-U))));
        ERRGE(3*(m-1)+n) = max(max(max(abs(UGE-U))));

    end
end
fprintf('                           Table 2\n');
fprintf('          Errors in solving the two-dimension model problem at t=1.\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('                        |              dt = 1/12          \n');
fprintf('         Method         | \t\t\th=0.1\t\t\t\tCPUTIME    \n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('   Backward-difference  | \t\t %.3e  \t\t %.3e    \n',ERRback(1),...
    ERRback(2));
fprintf('   Lawson-Morris        | \t\t %.3e  \t\t %.3e    \n',ERRLM(1),...
    ERRLM(2));
fprintf('   Global extrap        | \t\t %.3e  \t\t %.3e    \n',ERRGE(1),...
    ERRGE(2));
fprintf('-------------------------------------------------------------------------------------\n');