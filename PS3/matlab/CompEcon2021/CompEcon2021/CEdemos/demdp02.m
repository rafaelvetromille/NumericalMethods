%% DEMDP02 Asset Replacement Model
%
% Profit-maximizing manufacturer must decide when to replace an aging
% asset.

% States
%     p       profit contribution
%     a       asset age (1..A)
% Actions
%     j       keep(1) or replace(2) asset
% Parameters
%     A       maximum asset age 
%     alpha   production function coefficients
%     kappa   net replacement cost
%     pbar    long-run mean profit contribution
%     gamma   profit contribution autoregression coefficient
%     sigma   standard deviation of profit contribution shock
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';         % basis function type
n     =  100;               % number of collocation nodes
m     =    5;               % number of profit contribution shocks
nper  =   21;               % number of periods simulated
nrep  = 10000;              % number of replications

% Base case model parameters
A       = 6;                % maximum asset age 
alpha   = [50 -2.5 -2.5];   % production function coefficients
kappa   = 40;               % net replacement cost
pbar    = 1;                % long-run mean unit profit contribution
gamma   = 0.5;              % unit profit contribution autoregression coefficient
sigma   = 0.15;             % standard deviation of unit profit contribution shock
delta   = 0.9;              % discount factor 

% Profit contribution shock distribution
[e,w] = qnwnorm(m,0,sigma^2);

% Minimum and maximum profit contributions
pmin = pbar + min(e)/(1-gamma);
pmax = pbar + max(e)/(1-gamma);

% Approximation structure
basis = fundefn(basistype,n,pmin,pmax);     % basis functions

% Deterministic discrete state transitions
h = [2:A 1; ones(1,A)];

% Model structure
clear model
model.func = @func;                         % model functions
model.params = {A alpha kappa pbar gamma};  % function parameters
model.discount = delta;                   	% discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 0;                               % dimension of continuous action
model.ni = A;                               % number of discrete states
model.nj = 2;                               % number of discrete actions
model.e  = e;                              	% continuous state shocks
model.w  = w;                              	% continuous state shock probabilities
model.h  = h;                              	% deterministic discrete state transitions


%% BASE CASE MODEL SOLUTION

% Solve collocation equation
[p,v,~,pr,vr,~,resid] = dpsolve(model,basis); 

% Critical profit contributions and values
fprintf('Critical Replacement Profit Contributions\n')
pcrit = zeros(A,1);
vcrit = zeros(A,1);
for a=1:A-1
  pcrit(a) = interp1(vr(:,a,1)-vr(:,a,2),pr,0);
  vcrit(a) = interp1(pr,vr(:,a,1),pcrit(a));
  if isnan(pcrit), continue, end
  fprintf('   Age %2i  %5.2f\n'  ,a,pcrit(a))
end
fprintf('\n')

% Plot action-contingent value functions
figure
hold on
plot(pr,squeeze(vr(:,:,1)))
for a=1:A-1
  plotvdash(pcrit(a),vcrit(a))
  plotbullet(pcrit(a),vcrit(a))
  plottext(pcrit(a)+0.01,[],['$p^*_' int2str(a) '$'])
end
xtickformat('%.1f')
title('Action-Contingent Value Functions')
xlabel('Profit Contribution')
ylabel('Value of Production Operation')
legend('Keep a=1','Keep a=2','Keep a=3','Keep a=4','Keep a=5','Replace')

% Plot residual
figure
hold on
plot(pr,100*resid./max(vr,[],3))
plothdash([],0)
for a=1:A-1
  plotvdash(pcrit(a),0)
  plotbullet(pcrit(a),0)
  plottext(pcrit(a)+0.01,[],['$p^*_' int2str(a) '$'])
end
xtickformat('%.1f')
ytickformat('%.2f%%')
title('Bellman Equation Residual')
xlabel('Profit Contribution')
ylabel('Percent Residual')
legend('a=1','a=2','a=3','a=4','a=5','a=6','Location','NE')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm(0,sigma^2,nrep,nper);

% Initialize simulation
pinit = pbar*ones(nrep,1);              % initial profit contributions
ainit = ones(nrep,1);                   % initial asset ages

% Simulate model
[psim,~,asim,~] = dpsimul(model,basis,nper,pinit,ainit,pr,vr,[],esim);

% Ergodic moments
pavg = mean(psim(:));
aavg = mean(asim(:));
pstd = std(psim(:));
astd = std(asim(:));

% Print ergodic moments
fprintf('                     Ergodic      Ergodic\n') 
fprintf('                      Mean     Std Deviation\n') 
fprintf('Profit Contribution  %5.3f         %5.3f\n'  ,[pavg pstd])
fprintf('Age                  %5.3f         %5.3f\n\n',[aavg astd])

% Plot simulated and expected continuous state path
figure
hold on
plot(0:nper-1,psim(1:3,:))
plot(0:nper-1,mean(psim),'k')
ytickformat('%.1f')
title('Simulated and Expected Profit Contribution')
xlabel('Period')
ylabel('Profit Contribution')


%% PARAMETRIC PARAMETRIC SENSITIVITY ANALYSIS

% Intialization
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
aavgplot = zeros(nplot,1);

fprintf('Varying Net Replacement Cost\n')
modelplot = model;
kappamin = 30;
kappamax = 50;
kappaplot = nodeunif(nplot,kappamin,kappamax);
for ip=1:nplot
  rng('default')
  modelplot.params = {A alpha kappaplot(ip) pbar gamma};
  [~,~,~,sr,vr,~] = dpsolve(modelplot,basis);
  [psim,~,asim] = dpsimul(modelplot,basis,nper,pinit,ainit,sr,vr,[],esim);
  aavgplot(ip) = mean(asim(:,nper));
end
figure
plot(kappaplot,aavgplot)
ytickformat('%.1f')
xlabel('Net Replacement Cost')
ylabel('Mean Asset Age')

fprintf('Varying Long-Run Mean Profit Contribution\n')
modelplot = model;
pbarmin = 0.6;
pbarmax = 1.4;
pbarplot = nodeunif(nplot,pbarmin,pbarmax);
for ip=1:nplot
  rng('default')  
  pmin = pbarplot(ip) + min(e)/(1-gamma);
  pmax = pbarplot(ip) + max(e)/(1-gamma);
  basis = fundefn(basistype,n,pmin,pmax);
  modelplot.params = {A alpha kappa pbarplot(ip) gamma};
  [~,~,~,sr,vr,~] = dpsolve(modelplot,basis);
  [psim,~,asim] = dpsimul(modelplot,basis,nper,pinit,ainit,sr,vr,[],esim);
  aavgplot(ip) = mean(asim(:,nper));
end
figure
plot(pbarplot,aavgplot)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Long-Run Mean Profit Contribution')
ylabel('Mean Asset Age')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function out = func(flag,p,~,a,j,e,A,alpha,kappa,pbar,gamma)
switch flag
  case 'f'      % reward
    if j==2||a==A
      out = p*50-kappa;
    else
      out = p*(alpha(1)+alpha(2)*a+alpha(3)*a.^2);
    end
  case 'g'      % transition
    out = pbar + gamma*(p-pbar) + e;
end
end