% 2D plate in transient regime using explicit scheme
% Heat transfer - General program for any boundary conditions
% Created by: Jessica Guarato
% Last modified on: June, 2016

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
dt = 1;             % time step [s]
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

alfa = k/(rho*cp);
Fo = alfa*dt/(dx)^2;

T(:,:,1) = Ti*ones(ny,nx);
t(1) = 0;
nt = tfinal/dt+1;
balance(1) = 0;
Eac(1) = 0;

for p=2:nt
    
    t(p) = t(p-1)+dt;
    T(:,:,p) = T(:,:,p-1);
    
    for i=1:ny
        for j=1:nx
            
            % Calculation of temperatures
            % Top wall
            if i==1 & j>=2 & j<=nx-1
                if tp1==1
                    T(i,j,p) = T1;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*h1*dt/(rho*cp*dx))*(Fo*(2*T(i+1,j,p)+T(i,j-1,p)+T(i,j+1,p))+2*dt/(rho*cp*dx)*(h1*Tinf1+q1)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Bottom wall
            elseif i==ny & j>=2 & j<=nx-1
                if tp2==1
                    T(i,j,p) = T2;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*h2*dt/(rho*cp*dx))*(Fo*(2*T(i-1,j,p)+T(i,j-1,p)+T(i,j+1,p))+2*dt/(rho*cp*dx)*(h2*Tinf2+q2)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Left wall
            elseif  i>=2 & i<=ny-1 & j==1
                if tp3==1
                    T(i,j,p) = T3;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*h3*dt/(rho*cp*dx))*(Fo*(2*T(i,j+1,p)+T(i-1,j,p)+T(i+1,j,p))+2*dt/(rho*cp*dx)*(h3*Tinf3+q3)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Right wall
            elseif i>=2 & i<=ny-1 &  j==nx
                if tp4==1
                    T(i,j,p) = T4;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*h4*dt/(rho*cp*dx))*(Fo*(2*T(i,j-1,p)+T(i-1,j,p)+T(i+1,j,p))+2*dt/(rho*cp*dx)*(h4*Tinf4+q4)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Top left corner
            elseif i==1 & j==1
                if tp1==1
                    T(i,j,p) = T1;
                elseif tp3==1
                    T(i,j,p) = T3;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*dt/(rho*cp*dx)*(h1+h3))*(2*Fo*(T(i,j+1,p)+T(i+1,j,p))+2*dt/(rho*cp*dx)*(h1*Tinf1+h3*Tinf3+q1+q3)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Top right corner
            elseif i==1 & j==nx
                if tp1==1
                    T(i,j,p) = T1;
                elseif tp4==1
                    T(i,j,p) = T4;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*dt/(rho*cp*dx)*(h1+h4))*(2*Fo*(T(i,j-1,p)+T(i+1,j,p))+2*dt/(rho*cp*dx)*(h1*Tinf1+h4*Tinf4+q1+q4)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Bottom left corner
            elseif i==ny & j==1
                if tp2==1
                    T(i,j,p) = T2;
                elseif tp3==1
                    T(i,j,p) = T3;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*dt/(rho*cp*dx)*(h2+h3))*(2*Fo*(T(i,j+1,p)+T(i-1,j,p))+2*dt/(rho*cp*dx)*(h2*Tinf2+h3*Tinf3+q2+q3)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Bottom right corner
            elseif i==ny & j==nx
                if tp2==1
                    T(i,j,p) = T2;
                elseif tp4==1
                    T(i,j,p) = T4;
                else
                    T(i,j,p) = 1/(1+4*Fo+2*dt/(rho*cp*dx)*(h2+h4))*(2*Fo*(T(i,j-1,p)+T(i-1,j,p))+2*dt/(rho*cp*dx)*(h2*Tinf2+h4*Tinf4+q2+q4)+qv*dt/(rho*cp)+T(i,j,p-1));
                end
                
            % Rest of the plate
            else
                T(i,j,p) =1/(1+4*Fo)*(Fo*(T(i,j-1,p)+T(i,j+1,p)+T(i-1,j,p)+T(i+1,j,p))+qv*dt/(rho*cp)+T(i,j,p-1));
            end
        end
    end
end

% Energy balance

qconv1=0; qconv2=0; qconv3=0; qconv4=0; Eac(p)=0;
for i=1:ny
    for j=1:nx
        
        % Top wall
        if i==1
            if j>=2 & j<=nx-1
                qconv1 = h1*dx*(Tinf1-T(i,j,p))+qconv1;
                
            elseif j==1 || j==nx
                qconv1 = h1*dx/2*(Tinf1-T(i,j,p))+qconv1;
            end
        end
        
        % Bottom wall
        if i==ny
            if j>=2 & j<=nx-1
                qconv2 = h2*dx*(Tinf2-T(i,j,p))+qconv2;
                
            elseif j==1 || j==nx
                qconv2 = h2*dx/2*(Tinf2-T(i,j,p))+qconv2;
            end
        end
        
        % Left wall
        if j==1
            if i>=2 & i<=ny-1
                qconv3 = h3*dx*(Tinf3-T(i,j,p))+qconv3;
                
            elseif i==1 || i==ny & j==1
                qconv3 = h3*dx/2*(Tinf3-T(i,j,p))+qconv3;
            end
        end
        
        % Right wall
        if j==nx
            if i>=2 & i<=ny-1
                qconv4 = h4*dx*(Tinf4-T(i,j,p))+qconv4;
                
            elseif i==1 || i==ny
                qconv4 = h4*dx/2*(Tinf4-T(i,j,p))+qconv4;
            end
        end
        
        Eac(p) = Eac(p)+rho*cp*dx^2*(T(i,j,p)-T(i,j,p-1))/dt;
        
    end
end

balance(p) = qv*L*w+(q1+q2)*L+(q3+q4)*w+qconv1+qconv2+qconv3+qconv4-Eac(p);

% 2D graph of the temperature distribution on the plate
x = 0:dx:L;
y = 0:dy:w;

contourf(x,y,T(:,:,end))
set(gca,'ydir','reverse');
colorbar
xlabel('x (m)'); ylabel('y (m)'); title('Temperature Distribution [^oC]');