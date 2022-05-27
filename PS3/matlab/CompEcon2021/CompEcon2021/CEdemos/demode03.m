%% DEMODE03 Linear Initial Value Problem

%  Solve x1' =  1*x1 + 12*x2 - 60
%        x2' = -1*x1 -  6*x2 + 36
%  s.t.  x1(0)=5, x2(0)=2
%        t in [0,3]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'ODE Phase Diagram';
x1label  = '$x_1$';
x2label  = '$x_2$';

% Velocity Function
A  = [1 12; -1 -6];
b  = [-60; 36];
f = @(x) [A(1,1)*x(1,:)+A(1,2)*x(2,:)+b(1); ...
          A(2,1)*x(1,:)+A(2,2)*x(2,:)+b(2)];
        
% Closed-Form Solution 
X = @(t) [12 - 48*exp(-2*t) + 42*exp(-3*t) ...
           4 + 12*exp(-2*t) - 14*exp(-3*t)];
         
% Initial State
xinit = [6;2];
       
% Time Horizon
T = 3;
         
         
%% SOLVE ODE ANALYTICALLY
         
% Time Discretization
N = 100;              % number of time nodes
t = nodeunif(N,0,T);  % time nodes
         
% Plot Closed-Form Solution in Time Domain
figure
plot(t,X(t))
xtickformat('%.1f')
title('Solution')
xlabel('Time')
legend(x1label,x2label)


%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Solve ODE
[x,t] = oderk4(f,xinit,T,N);

% Plot Runge-Kutta Approximation Errors
figure
plot(t,x-X(t),t,0*t,'k:')
xtickformat('%.1f')
title('Runge Kutta Approximation Errors')
xlabel('Time')
ylabel('Error')
legend(x1label,x2label)


%% SOLVE ODE USING COLLOCATION METHOD (ODECOL)

% Solve ODE
n  = 20;    % number of basis functions
[x,t,res] = odecol(f,xinit,T,n);

% Plot Collocation Approximation Errors
figure
plot(t,x-X(t),t,0*t,'k:')
xtickformat('%.1f')
ytickformat('%.1f')
title('Collocation Approximation Errors')
xlabel('Time')
ylabel('Error')
legend(x1label,x2label)

% Plot Residuals
figure
hold on
plot(t,res)
plothdash([],0)
xtickformat('%.1f')
ytickformat('%.1f')
title('Collocation Residuals')
xlabel('Time')
ylabel('Residual')
legend(x1label,x2label)


%% STEADY-STATE

% Compute Steady State
xstst = -A\b;
disp('Steady State')
disp(xstst)
disp('Eigenvalues')
disp(eig(fdjac(f,xstst)))


%% PHASE DIAGRAM

% Plotting Limits
x1lim = [0 15];  % x1 plotting limits
x2lim = [0  8];  % x2 plotting limits

% Compute Separatrix
xspx = odespx(f,xstst,12);

% Compute Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = [x1 -(A(1,1)*x1+b(1))/A(1,2)]; 
x2null = [x1 -(A(2,1)*x1+b(2))/A(2,2)]; 

% Plot Phase Diagram
odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,xspx,xstst,x1null,x2null)


%% SAVE FIGURES
printfigures(mfilename)