% 2D plate in transient regime using implicit scheme
% Heat transfer - General program for any boundary conditions
% Made by: Jessica Guarato
% Last modified on: March 08, 2018

clc; clear all; close all;

%% Setup

% Problem conditions
L = 0.5;            % x [m]
w = 0.25;           % y [m]
k = 15;             % thermal conductivity [W/(m.k)]
rho = 7900;         % density [kg/m3]
cp = 477;           % specific heat[J/(kg.K)]
qv = 1e4;           % volumetric heat generation [W/m^3]
nx = 31;            % number of nodes in x (must be an odd number)
dt = 100;             % time step [s]
tfinal = 1e6;       % final time step [s]

% Initial condition
Ti = 20;

% Top wall (1)
tp1 = 0;            % [1] if prescribed temperature on the wall; [0] if not
T1 = 0;             % prescribed temperature on the wall [degrees Celsius]
h1 = 0;             % convection coefficient, if it has convective fluid in contact with the wall
Tinf1 = 0;          % convective fluid temperature [degrees Celsius]
q1 = 1e3;           % heat flow on the wall [W/m^2]

% Bottom wall (2)
tp2 = 0;            % [1] if prescribed temperature on the wall; [0] if not
T2 = 0;             % prescribed temperature on the wall [degrees Celsius]
h2 = 0;             % convection coefficient, if it has convective fluid in contact with the wall
Tinf2 = 0;          % convective fluid temperature [degrees Celsius]
q2 = 0;             % heat flow on the wall [W/m^2]

% Left wall (3)
tp3 = 0;            % [1] if prescribed temperature on the wall; [0] if not
T3 = 0;             % prescribed temperature on the wall [degrees Celsius]
h3 = 50;            % convection coefficient, if it has convective fluid in contact with the wall
Tinf3 = 10;         % convective fluid temperature [degrees Celsius]
q3 = 0;             % heat flow on the wall [W/m^2]

% Right wall (4)
tp4 = 0;            % [1] if prescribed temperature on the wall; [0] if not
T4 = 0;             % prescribed temperature on the wall [degrees Celsius]
h4 = 50;            % convection coefficient, if it has convective fluid in contact with the wall
Tinf4 = 10;         % convective fluid temperature [degrees Celsius]
q4 = 0;             % heat flow on the wall [W/m^2]

%% Problem solution

% Obtaining parameters
dx = L/(nx-1);
dy = dx;
ny = w/dx+1;        % number of nodes in y
x = linspace(0,L,nx);
y = linspace(0,w,ny);

alfa = k/(rho*cp);
Fo = alfa*dt/(dx)^2;

T = (Ti*ones(nx*ny))';
Tt(:,:,1) = Ti*ones(ny,nx);
t(1) = 0;
nt = tfinal/dt+1;

for p=2:nt
    
    A = zeros(nx*ny,nx*ny);
    
    % Top right corner
    for z=1
        A(z,z) = -(4*Fo+2*dx*Fo/k*h1+2*h3*dx/k*Fo+1);
        A(z,z+1) = 2*Fo;
        A(z,z+nx) = 2*Fo;
        B(z) = -(2*dx*Fo/k)*(q1+q3+h1*Tinf1+h3*Tinf3)-qv*(dt/rho/cp)-T(z);
    end
    
    % Top wall
    for z=2:nx-1
        A(z,z-1) = Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h1+1);
        A(z,z+1) = Fo;
        A(z,z+nx) = 2*Fo;
        B(z) = -(2*dx*Fo/k)*(q1+h1*Tinf1)-qv*(dt/rho/cp)-T(z);
    end
    
    % Top right corner
    for z=nx
        A(z,z-1) = 2*Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h1+2*h4*dx/k*Fo+1);
        A(z,z+nx) = 2*Fo;
        B(z) = -(2*dx*Fo/k)*(q1+q4+h1*Tinf1+h4*Tinf4)-qv*(dt/rho/cp)-T(z);
    end
    
    % Left wall
    for z=nx+1:nx:(nx*ny-nx-1)
        A(z,z-nx) = Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h3+1);
        A(z,z+1) = 2*Fo;
        A(z,z+nx) = Fo;
        B(z) = -(2*dx*Fo/k)*(q3+h3*Tinf3)-qv*(dt/rho/cp)-T(z);
    end
    
    % Interior
    for f=2:ny-1
        for z=(f*nx-nx+2):(f*nx-1)
            A(z,z-nx) = Fo;
            A(z,z-1)= Fo;
            A(z,z) = -(4*Fo+1);
            A(z,z+1) = Fo;
            A(z,z+nx) = Fo;
            B(z) = -qv*(dt/rho/cp)-T(z);
        end
    end
    
    % Right wall
    for z=2*nx:nx:(nx*ny-nx)
        A(z,z-nx) = Fo;
        A(z,z-1) = 2*Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h4+1);
        A(z,z+nx) = Fo;
        B(z) = -(2*dx*Fo/k)*(q4+h4*Tinf4)-qv*(dt/rho/cp)-T(z);
    end
    
    % Bottom left corner
    for z=nx*ny-nx+1
        A(z,z-nx) = 2*Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h2+2*h3*dx/k*Fo+1);
        A(z,z+1) = 2*Fo;
        B(z) = -(2*dx*Fo/k)*(q2+q3+h2*Tinf2+h3*Tinf3)-qv*(dt/rho/cp)-T(z);
    end
    
    % Bottom wall
    for z=(nx*ny-nx+2):(nx*ny-1)
        A(z,z-nx) = 2*Fo;
        A(z,z-1) = Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h2+1);
        A(z,z+1) = Fo;
        B(z) = -(2*dx*Fo/k)*(q2+h2*Tinf2)-qv*(dt/rho/cp)-T(z);
    end
    
    % Bottom right corner
    for z=nx*ny
        A(z,z-nx) = 2*Fo;
        A(z,z-1) = 2*Fo;
        A(z,z) = -(4*Fo+2*dx*Fo/k*h2+2*h4*dx/k*Fo+1);
        B(z) = -(2*dx*Fo/k)*(q2+q4+h2*Tinf2+h4*Tinf4)-qv*(dt/rho/cp)-T(z);
    end
    
    C = B';
    T = A\C;
    
    for ff=1:ny
        for z=(ff*nx-nx+1):(ff*nx)
            fff = z-(ff-1)*nx;
            Tt(ff,fff,p) = T(z);
        end
    end
    
    % 2D graph of the temperature distribution on the plate
    contourf(x,y,Tt(:,:,p))
    set(gca,'ydir','reverse');
    colorbar
    xlabel('x [m]'); ylabel('y [m]'); title('Temperature Distribution [^oC]');
    pause(0.1);
    
end

%% Energy balance

balance(1) = 0;
Eac(1) = 0;

for p=1:nt-1
    
    qconv1=0; qconv2=0; qconv3=0; qconv4=0;
    Eac(p+1)=0;
    
    for i=1:ny
        for j=1:nx
            
            % Interior
            eac = rho*cp*dx^2*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
            
            % Top wall
            if i==1
                if j>=2 & j<=nx-1
                    qconv1 = h1*dx*(Tinf1-Tt(i,j,p+1))+qconv1;
                    eac = rho*cp*(dx^2/2)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                    
                elseif j==1 || j==nx
                    qconv1 = h1*dx/2*(Tinf1-Tt(i,j,p+1))+qconv1;
                    eac = rho*cp*(dx^2/4)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                end
            end
            
            % Bottom wall
            if i==ny
                if j>=2 & j<=nx-1
                    qconv2 = h2*dx*(Tinf2-Tt(i,j,p+1))+qconv2;
                    eac = rho*cp*(dx^2/2)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                    
                elseif j==1 || j==nx
                    qconv2 = h2*dx/2*(Tinf2-Tt(i,j,p+1))+qconv2;
                    eac = rho*cp*(dx^2/4)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                end
            end
            
            % Left wall
            if j==1
                if i>=2 & i<=ny-1
                    qconv3 = h3*dx*(Tinf3-Tt(i,j,p+1))+qconv3;
                    eac = rho*cp*(dx^2/2)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                    
                elseif i==1 || i==ny
                    qconv3 = h3*dx/2*(Tinf3-Tt(i,j,p+1))+qconv3;
                end
            end
            
            % Right wall
            if j==nx
                if i>=2 & i<=ny-1
                    qconv4 = h4*dx*(Tinf4-Tt(i,j,p+1))+qconv4;
                    eac = rho*cp*(dx^2/2)*(Tt(i,j,p+1)-Tt(i,j,p))/dt;
                    
                elseif i==1 || i==ny
                    qconv4 = h4*dx/2*(Tinf4-Tt(i,j,p+1))+qconv4;
                end
            end
            
            Eac(p+1) = Eac(p+1)+eac;
            
        end
    end
    
    balance(p+1) = qv*L*w+(q1+q2)*L+(q3+q4)*w+qconv1+qconv2+qconv3+qconv4-Eac(p+1);
    
end

% To check if the largest balance value (as a function of time) is close to zero:
max_balance = max(abs(balance)) 