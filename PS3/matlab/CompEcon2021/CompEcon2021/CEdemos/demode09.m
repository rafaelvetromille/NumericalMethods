%% DEMODE09 Tobin's Q Model

%  Solve
%    k' = k.*(max(q-1,0)/b-delta)
%    q' = (delta+rho)*q-1./k-(q-1)^2/(2*b)]
%  where
%    k: capital stock
%    q: price of capital
%    x = [k;q]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'Tobin''s Q Model Phase Diagram';
x1label  = 'Capital Stock';
x2label  = 'Price of Capital';

% Velocity Function
delta = 0.05;      % depreciation rate
b     = 5;         % marginal adjustment cost
rho   = 0.05;      % interest rate
k = @(x) x(1,:);
q = @(x) x(2,:);
f = @(x) [k(x).*(max(q(x)-1,0)/b-delta); ...
  (delta+rho)*q(x)-1./k(x)-(q(x)-1).^2/(2*b)];

% Initial States
xinit = [9 9 9; 1.1 1.2074 1.3];
       
% Time Horizons
T = [10 30 10];


%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Solve for Different Initial Values
N = 1000;  % number of time nodes
x = zeros(N,2,3);
for i=1:3
  [x(:,:,i),t] = oderk4(f,xinit(:,i),T(i),N);
end

% Plot Solutions in Time Domain
for i=1:3
  figure
  plot(t,x(:,:,i))
  legend(x1label,x2label)
end


%% STEADY-STATE

% Compute Steady State
qstst = 1+b*delta;
kstst = 1/((delta+rho)*qstst-(b*delta)^2/(2*b));
xstst = [kstst;qstst];
disp('Steady State')
disp(xstst)
disp('Eigenvalues')
disp(eig(fdjac(f,xstst)))


%% PHASE DIAGRAM

% Ploting Limits
x1lim = [  7  10];  % k plotting limit
x2lim = [1.0 1.5];  % q plotting limit

% Compute Separatrix
xspx = odespx(f,xstst,300);

% Compute Nullclines
k = nodeunif(100,x1lim(1),x1lim(2));
A = -1/(2*b);
B = (delta+rho)+1/b;
C = 1./k+1/(2*b);
qnull = [k real(-B+sqrt(B^2+4*A*C))/(2*A)];
knull = [[kstst;kstst] x2lim'];

% Plot Phase Diagram
odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,xspx,xstst,qnull,knull)
xtickformat('%.1f')
ytickformat('%.1f')

% Label State Paths
for i=1:3
  text(x(1,1,i),x(1,2,i),['\bf' num2str(i)])
end


%% SAVE FIGURES
printfigures(mfilename)