%% DEMODE06 Predator-Prey Model

%  Solve
%    x1' = -a1*x1         + a2*x1*x2 - h
%    x2' =        + b1*x2 - b2*x1*x2
%  where
%    x1: Predator population
%    x2: Prey population

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'Predator-Prey Model Phase Diagram';
x1label  = 'Predator';
x2label  = 'Prey';

% Model Parameters
a1 = 0.40;      % natural death rate of predators
a2 = 0.01;      % growth rate of predators per predation
b1 = 1.00;      % natural growth rate of prey
b2 = 0.02;      % death rate of prey per predation

% Initial State
xinit = [40;40];
       
% Time Horizon
T = 30;


%% PREDATOR-PREY MODEL WITHOUT (H=0) AND WITH (H=3) HUNTING

for h=[0 3]
  
  % Velocity Function
  f = @(x) [-a1*x(1,:)+a2*x(1,:).*x(2,:)-h;b1*x(2,:)-b2*x(1,:).*x(2,:)];
  
  % Solve ODE Using Collocation Method
  n = 101;          % number of basis functions
  [x,t] = odecol(f,xinit,T,n);
  
  % Plot Solution in Time Domain
  figure
  plot(t,x)
  title('Population Paths')
  xlabel('Time')
  legend(x1label,x2label)
  
  % Compute Steady State
  xstst = [b1/b2;(h*b2+a1*b1)/(a2*b1)];
  disp('Steady State')
  disp(xstst)
  disp('Eigenvalues')
  disp(eig(fdjac(f,xstst)))
  
  % Plotting Limits
  x1lim = [0  80];         % Predator plotting limits
  x2lim = [0 100];         % Prey plotting limits
  
  % Compute Nullclines
  x1 = nodeunif(100,x1lim(1),x1lim(2));
  x1null = [x1 (a1*x1+h)./(a2*x1)];
  x2 = nodeunif(100,x2lim(1),x2lim(2));
  x2null = [ones(100,1)*b1/b2 x2];
  
  % Plot Phase Diagram
  odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,[],xstst,x1null,x2null)
  
end


%% SAVE FIGURES
printfigures(mfilename)