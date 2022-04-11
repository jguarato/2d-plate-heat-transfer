% 2D plate in steady state
% Heat transfer - General program for any boundary conditions
% Made by: Jessica Guarato
% Last modified on: June, 2016

clc; close all; clear all;

%% Setup

L = 0.5;            % x [m]
w = 0.25;           % y [m]
k = 15;             % thermal conductivity [W/(m.k)]
tolerance = 1e-8;   % tolerance of (current iteration - previous iteration)
nx = 31;            % number of nodes in x (must be an odd number)
qv = 1e4;           % volumetric heat generation [W/m^3]

% Top wall(1)
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

T = zeros(ny,nx);
max_error = 1;

while max_error>tolerance
    
    To = T;
    
    for i=1:ny
        for j=1:nx
            
            % Calculation of temperatures
            
            % Top wall
            if i==1 & j>=2 & j<=nx-1
                if tp1==1
                    T(i,j) = T1;
                else
                    T(i,j) = 1/(2*k/dx+h1)*(k/(2*dx)*(T(i,j-1)+T(i,j+1))+k/dx*T(i+1,j)+q1+h1*Tinf1+dx/2*qv);
                end
                
            % Bottom wall
            elseif i==ny & j>=2 & j<=nx-1
                if tp2==1
                    T(i,j) = T2;
                else
                    T(i,j) = 1/(2*k/dx+h2)*(k/(2*dx)*(T(i,j-1)+T(i,j+1))+k/dx*T(i-1,j)+q2+h2*Tinf2+dx/2*qv);
                end
                
            % Left wall
            elseif  i>=2 & i<=ny-1 & j==1
                if tp3==1
                    T(i,j) = T3;
                else
                    T(i,j) = 1/(2*k/dx+h3)*(k/(2*dx)*(T(i-1,j)+T(i+1,j))+k/dx*T(i,j+1)+q3+h3*Tinf3+dx/2*qv);
                end
                
            % Right wall
            elseif i>=2 & i<=ny-1 &  j==nx
                if tp4==1
                    T(i,j) = T4;
                else
                    T(i,j) = 1/(2*k/dx+h4)*(k/(2*dx)*(T(i-1,j)+T(i+1,j))+k/dx*T(i,j-1)+q4+h4*Tinf4+dx/2*qv);
                end
                
            % Top left corner
            elseif i==1 & j==1
                if tp1==1
                    T(i,j) = T1;
                elseif tp3==1
                    T(i,j) = T3;
                else
                    T(i,j) = 1/(2*k/dx+h1+h3)*(k/dx*(T(i,j+1)+T(i+1,j))+h1*Tinf1+h3*Tinf3+q1+q3+dx/2*qv);
                end
                
            % Top right corner
            elseif i==1 & j==nx
                if tp1==1
                    T(i,j) = T1;
                elseif tp4==1
                    T(i,j) = T4;
                else
                    T(i,j) = 1/(2*k/dx+h1+h4)*(k/dx*(T(i,j-1)+T(i+1,j))+h1*Tinf1+h4*Tinf4+q1+q4+dx/2*qv);
                end
                
            % Bottom left corner
            elseif i==ny & j==1
                if tp2==1
                    T(i,j) = T2;
                elseif tp3==1
                    T(i,j) = T3;
                else
                    T(i,j) = 1/(2*k/dx+h2+h3)*(k/dx*(T(i,j+1)+T(i-1,j))+h2*Tinf2+h3*Tinf3+q2+q3+dx/2*qv);
                end
                
            % Bottom right corner
            elseif i==ny & j==nx
                if tp2==1
                    T(i,j) = T2;
                elseif tp4==1
                    T(i,j) = T4;
                else
                    T(i,j) = 1/(2*k/dx+h2+h4)*(k/dx*(T(i,j-1)+T(i-1,j))+h2*Tinf2+h4*Tinf4+q2+q4+dx/2*qv);
                end
                
            % Rest of the plate
            else
                T(i,j) = (T(i,j+1)+T(i,j-1)+T(i+1,j)+T(i-1,j)+dx^2/k*qv)/4;
            end
            
        end
    end
    
    error = abs(T-To);
    max_error = max(max(error));
end

% 2D graph of the temperature distribution on the plate
x = linspace(0,L,nx);
y = linspace(0,w,ny);

contourf(x,y,T)
set(gca,'ydir','reverse');
colorbar
xlabel('x [m]'); ylabel('y [m]'); title('Temperature Distribution [^oC]');

%% Energy balance

qconv1=0; qconv2=0; qconv3=0; qconv4=0;
for i=1:ny
    for j=1:nx
        
        % Top wall
        if i==1
            if j>=2 & j<=nx-1
                qconv1 = h1*dx*(Tinf1-T(i,j))+qconv1;
                
            elseif j==1 || j==nx
                qconv1 = h1*dx/2*(Tinf1-T(i,j))+qconv1;
            end
        end
        
        % Bottom wall
        if i==ny
            if j>=2 & j<=nx-1
                qconv2 = h2*dx*(Tinf2-T(i,j))+qconv2;
                
            elseif j==1 || j==nx
                qconv2 = h2*dx/2*(Tinf2-T(i,j))+qconv2;
            end
        end
        
        % Left wall
        if j==1
            if i>=2 & i<=ny-1
                qconv3 = h3*dx*(Tinf3-T(i,j))+qconv3;
                
            elseif i==1 || i==ny & j==1
                qconv3 = h3*dx/2*(Tinf3-T(i,j))+qconv3;
            end
        end
        
        % Right wall
        if j==nx
            if i>=2 & i<=ny-1
                qconv4 = h4*dx*(Tinf4-T(i,j))+qconv4;
                
            elseif i==1 || i==ny
                qconv4 = h4*dx/2*(Tinf4-T(i,j))+qconv4;
            end
        end
    end
end

balance = qv*L*w+(q1+q2)*L+(q3+q4)*w+qconv1+qconv2+qconv3+qconv4 