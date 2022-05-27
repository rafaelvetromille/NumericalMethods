%% DEMODE04 - Non-IVP Non-Homogeneous Linear ODE Example

%  Solve x1' = -1*x1 - 0.5*x2 + 2
%        x2' =        -0.5*x2 + 1
%  s.t   x1(0)=1, x2(1)=1
%        t in [0,10]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'ODE Phase Diagram';
x1label  = '$x_1$';
x2label  = '$x_2$';

% Velocity Function
A = [-1 -0.5; 0 -0.5];
b = [2; 1];
f = @(x) [A(1,1)*x(1,:)+A(1,2)*x(2,:)+b(1); ...
          A(2,1)*x(1,:)+A(2,2)*x(2,:)+b(2)];
        
% Boundary Conditions
bx = [1;2];     % boundary variables
bt = [0;1];     % boundary times
bv = [1;1];     % boundary values

% Time Horizon
T = 10;

% Closed-Form Solution 
X = @(t) [1-exp(1/2-t)+exp((1-t)/2) 2-exp((1-t)/2)];


%% SOLVE ODE ANALYTICALLY

% Time Discretization
N = 200;              % number of time nodes
t = nodeunif(N,0,T);  % time nodes

% Plot Closed-Form Solution in Time Domain
figure
plot(t,X(t))
ytickformat('%.1f')
title('Solution')
xlabel('Time')
legend(x1label,x2label)


%% SOLVE ODE USING COLLOCATION METHOD (ODECOL)

% Solve ODE
n = 15;    % number of basis functions
[x,t,res] = odecol(f,bv,T,n,bt,bx);

% Plot Collocation Approximation Errors
figure
plot(t,x-X(t),t,0*t,'k:')
title('Collocation Approximation Errors')
xlabel('Time')
ylabel('Error')
legend(x1label,x2label)

% % Plot Residuals
% figure
% hold on
% plot(t,res)
% plothdash([],0)
% title('Collocation Residuals')
% xlabel('Time')
% ylabel('Residual')
% legend(x1label,x2label)


%% STEADY-STATE

% Compute Steady State
xstst = -A\b;
disp('Steady State')
disp(xstst)
disp('Eigenvalues')
disp(eig(fdjac(f,xstst)))


%% PHASE DIAGRAM

% Plotting Limits
x1lim = [0 2];  % x1 plotting limits
x2lim = [0 4];  % x2 plotting limits

% Compute Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = [x1 -(A(1,1)*x1+b(1))/A(1,2)]; 
x2null = [x1 -(A(2,1)*x1+b(2))/A(2,2)]; 

% Plot Phase Diagram
odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,[],xstst,x1null,x2null)
xtickformat('%.1f')


%% SAVE FIGURES
printfigures(mfilename)