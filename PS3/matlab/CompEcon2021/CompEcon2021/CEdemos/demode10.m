%% DEMODE10 Regional Migration Model

%  Solve
%    L' = L*sqrt(max(0,B-k0)/alpha)
%    B' = rho*B-L.^(-theta)+wbar
%  where
%    L: labor
%    B: benefit from migration
%    x = [L;B]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'Regional Migration Model Phase Diagram';
x1label  = 'Labor';
x2label  = 'Benefit';

% Velocity Function
theta = 0.5;          % inverse labor demand elasticity
k0    = 0.1;          % minimum migration cost
k1    = 0.1;          % migration cost parameter
k2    = 0.1;          % migration cost parameter
rho   = 0.05;         % discount rate
wbar  = 1;            % world wage rate
L = @(x) x(1,:);
B = @(x) x(2,:);
f = @(x) [L(x).*(sqrt(k1^2+4*k2*(B(x)-k0))-k1)/(2*k2); rho*B(x)-L(x).^(-theta)+wbar];

% Initial States
xinit = [ 0.8  0.8  0.8; 0.12 0.1542 0.19];
       
% Time Horizon
T = 2;


%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Solve for Different Initial Values
N = 1000;   % number of time nodes
x = zeros(N,2,3);
for i=1:3
  [x(:,:,i),t] = oderk4(f,xinit(:,i),T,N);
end

% Plot Solutions in Time Domain
for i=1:3
  figure
  plot(t,x(:,:,i))
  xtickformat('%.1f')
  ytickformat('%.1f')
  legend(x1label,x2label)
end


%% STEADY-STATE

% Compute Steady State
Bst = k0;
Lst = (wbar+rho*Bst)^(-1/theta);
xstst = [Lst;Bst];
disp('Steady State')
disp(xstst)
disp('Eigenvalues')
disp(eig(fdjac(f,xstst)))


%% PHASE DIAGRAM

% Ploting Limits
Llim = [0.60 1.40];   % labor plotting limits
Blim = [0.05 0.20];   % benefit plotting limits

% Compute Separatrix
xspx = odespx(f,xstst,100);
i = find(xspx(:,2)>Bst);
xspx = xspx(i,:);

% Compute Nullclines
L = nodeunif(100,Llim(1),Llim(2));
Bnull = [L (L.^(-theta)-wbar)/rho];
Lnull = [Llim' [Bst;Bst]];

% Plot Phase Diagram
odephase(Llim,Blim,figtitle,x1label,x2label,f,x,xspx,xstst,Bnull,Lnull)
xtickformat('%.1f')
ytickformat('%.2f')

% Label State Paths
for i=1:3
  text(x(1,1,i),x(1,2,i),['\bf' num2str(i)])
end


%% SAVE FIGURES
printfigures(mfilename)