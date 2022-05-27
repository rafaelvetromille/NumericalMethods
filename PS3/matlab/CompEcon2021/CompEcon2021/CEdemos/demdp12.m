%% DEMDP12 Production Management Model
%
% Profit maximizing entrepeneur must decide how much to produce, subject to 
% production adjustment costs.

% States
%     i       market price (discrete)
%     s       lagged production (continuous)
% Actions
%     x       current production
% Parameters
%     alpha   marginal adjustment cost
%     beta    marginal production cost parameters
%     pbar    long-run average market price
%     sigma   market price volatility
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)



%% FORMULATION

% Numerical control parameters
basistype = 'cheb';                         % basis function type
n =      21;                                % degree of approximation
smin =    0;                                % minimum production
smax =   20;                                % maximum production
m    =    3;                                % number of shocks
nper =   31;                                % number of periods simulated
nrep = 2000;                                % number of replications
nplot =  31;                                % number of parameter values plotted

% Base case model parameters
alpha = 0.01;                               % marginal adjustment cost
beta  = [0.8 0.03];                         % marginal production cost parameters
pbar  = 1.0;                                % long-run average market price
sigma = 0.2;                                % market price shock standard deviation
delta = 0.9;                                % discount factor

% Approximation structure
basis = fundefn(basistype,n,smin,smax);

% Continuous state shock distribution
[p,w] = qnwlogn(m,pbar,sigma^2);            % shocks and probabilities
q = w(:,ones(1,m))';                        % transition probabilities

% Model structure
clear model
model.func = @func;                         % model functions
model.params = {alpha beta p};              % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 1;                               % dimension of continuous action
model.ni = m;                               % number of discrete states
model.nj = 0;                               % number of discrete actions
model.q  = q;                               % discrete state transition probabilities


%% BASE CASE MODEL SOLUTION

% Nonstochatic steady-state
sstar = (pbar-beta(1))/beta(2);

% Check model derivatives
dpcheck(model,sstar,sstar)

% Solve collocation equation
[v,x,c,sr,vr,xr,resid] = dpsolve(model,basis);

% Plot optimal policy
figure
plot(sr,xr)
title('Optimal Production Policy')
xlabel('Lagged Production')
ylabel('Production')
legend('Low Price','Average Price','High Price')

% Plot value function
figure
plot(sr,vr)
title('Value Function')
xlabel('Lagged Production')
ylabel('Value of the Firm')
legend('Low Price','Average Price','High Price')

% Plot shadow price function
figure
lambda = funeval(c,basis,sr,1);
plot(sr,lambda)
ytickformat('%.2f')
title('Shadow Price of Lagged Production')
xlabel('Lagged Production')
ylabel('Shadow Price')
legend('Low Price','Average Price','High Price')

% Plot residual
figure
hold on
plot(sr,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Lagged Production')
ylabel('Residual')
legend('Low Price','Average Price','High Price')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(pbar,sigma^2,nrep,nper);

% Initialize simulation
sinit = smin(ones(nrep,1),:);           % initial lagged production
iinit = 2*ones(nrep,1);                 % initial market price state

% Simulate model
[ssim,xsim,isim] = dpsimul(model,basis,nper,sinit,iinit,sr,vr,xr,esim);
psim = p(isim);

% Ergodic moments
savg = mean(ssim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
pstd = std(psim(:)); 

% Print ergodic moments
fprintf('\n')
fprintf('Ergodic Moments\n') 
fprintf('          Nonstochastic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Price        %5.3f         %5.3f         %5.3f\n'  ,[pbar  pavg pstd])
fprintf('Production   %5.3f         %5.3f         %5.3f\n\n',[sstar savg sstd])

% Plot simulated action path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Production')
xlabel('Period')
ylabel('Production')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Intialization
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
xavgplot = zeros(nplot,1);
xstdplot = zeros(nplot,1);
vinit = v;
xinit = x;

fprintf('Varying Marginal Adjustment Cost\n')
modelplot = model;
alphamin = 0.0;
alphamax = 0.2;
alphaplot = nodeunif(nplot,alphamin,alphamax);
v = vinit; x = xinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {alphaplot(ip) beta p};
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,v,x);
  [ssim,xsim,isim] = dpsimul(modelplot,basis,nper,sinit,iinit,sr,vr,xr,esim);
  xavgplot(ip) = mean(xsim(:,nper));
  xstdplot(ip) = std(xsim(:,nper));
end
figure
plot(alphaplot,xavgplot)
xtickformat('%.2f')
ytickformat('%.1f')
xlabel('Marginal Adjustment Cost')
ylabel('Ergodic Mean Production')
figure
plot(alphaplot,xstdplot)
xtickformat('%.2f')
xlabel('Marginal Adjustment Cost')
ylabel('Ergodic Standard Deviation of Production')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,q,i,~,~,alpha,beta,p)
n = length(s);
p = p(i);
l = s;
switch flag
  case 'b'      % bounds
    out1 = zeros(n,1);
    out2 = inf*ones(n,1);
    out3 = [];
  case 'f'      % reward
    out1 = p.*q - (beta(1)*q+0.5*beta(2)*q.^2) - 0.5*alpha*((q-l).^2);
    out2 = p - beta(1) - beta(2)*q - alpha*(q-l);
    out3 = (-beta(2)-alpha)*ones(n,1);
  case 'g'      % transition
    out1 = q;
    out2 = ones(n,1);
    out3 = zeros(n,1);
end
end