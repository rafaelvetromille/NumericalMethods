%% DEMREM01 Asset Pricing Model
%
% Lucas-Prescott rational expectations asset pricing model.

% States
%     d       asset dividend
% Response
%     p       asset price
% Parameters
%     alpha    coefficient of risk aversion
%     dbar    long-run mean dividend
%     gamma   dividend autoregression coefficient
%     sigma   dividend shock standard deviation
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';                     % basis function type
n      = 25;                            % degree of approximation
dmin   = 0.1;                           % minimum state
dmax   = 1.9;                           % maximum state
m      = 5;                             % number of shocks
nper   = 50;                            % number of periods simulated
nrep   = 1000;                          % number of replications
nplot = 35;                             % number of parameter values plotted

% Base case model parameters
alpha = 0.5;                            % coefficient of risk aversion
dbar  = 1.0;                            % long-run mean dividend
gamma = 0.5;                            % dividend autoregression coefficient
sigma = 0.1;                            % dividend shock standard deviation
delta = 0.9;                            % discount factor

% Continuous state shock distribution
[e,w] = qnwnorm(m,0,sigma^2);           % shocks and weights

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {delta alpha dbar gamma};% function parameters
model.e  = e;                           % shocks
model.w  = w;                           % shock probabilities

% Approximation structure
[basis,Phi,dnode]  = fundefn(basistype,n,dmin,dmax);   % basis functions


%% BASE CASE MODEL SOLUTION

% Nonstochastic steady-state
pstar = delta*dbar/(1-delta);                       % asset price

% Solve rational expectations equilibrium directly
tic
LHS = diag(dnode.^(-alpha))*Phi;
RHS = 0;
for k=1:m
  dnext = dbar + gamma*(dnode-dbar) + e(k);
  LHS   = LHS - delta*w(k)*diag(dnext.^(-alpha))*funbase(basis,dnext);
  RHS   = RHS + delta*w(k)*dnext.^(1-alpha);
end
c = LHS\RHS;
p = funeval(c,basis,dnode);

% Residual
dd = nodeunif(10*n,dmin,dmax);
Ef = 0;
for k=1:m
  dnext = dbar + gamma*(dd-dbar) + e(k);
  f     = diag(dnext.^(-alpha))*(funeval(c,basis,dnext)+dnext);
  Ef    = Ef + delta*w(k)*f;
end
toc

% Solve rational expectations equilibrium using remsolve
[p,c,dr,pr,fr,resid] = remsolve(model,basis,p);

% Plot equilibrium asset price
figure
plot(dr,pr)
xtickformat('%.1f')
title('Equilibrium Asset Price')
xlabel('Dividend')
ylabel('Price')  

% Plot expected arbitrage profit function
figure
hold on
plot(dr,fr)
plothdash([],0)
xtickformat('%.1f')
title('Expected Arbitrage Benefit')
xlabel('Dividend')
ylabel('Profit')

% Plot residual
figure
hold on
plot(dr,resid)
plothdash([],0)
xtickformat('%.1f')
ytickformat('%.1f')
title('Response Residual')
xlabel('Dividend')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm(0,sigma^2,nrep,nper+1);

% Initialize simulation
dinit = ones(nrep,1);
rng('default')

% Preallocate arrays
dsim = zeros(nrep,nper);
psim = zeros(nrep,nper);

% Simulate model directly, without remsimul
dd = dinit;
for ip=1:nper+1
  pp = interp1(dr,pr,dd,'PCHIP');
  dsim(:,ip,:) = dd;
  psim(:,ip,:) = pp;
  if ip<nper+1
    dd = dbar+gamma*(dd-dbar)+e(discrand(nrep,w),:);
  end
end

% Simulate model using remsimul
rng('default')
[dsim,psim] = remsimul(model,basis,nper,dinit,dr,pr,esim);

% Plot simulated and expected state path
figure
plot(0:nper,dsim(1:3,:),0:nper,mean(dsim),'k')
ytickformat('%.1f')
title('Simulated and Expected Asset Dividend')
xlabel('Period')
ylabel('Dividend')

% Plot simulated and expected response path
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
ytickformat('%.1f')
title('Simulated and Expected Asset Price')
xlabel('Period')
ylabel('Price')

% Ergodic moments
davg = mean(dsim(:)); 
pavg = mean(psim(:)); 
dstd = std(dsim(:)); 
pstd = std(psim(:)); 

% Print ergodic moments
fprintf('\n') 
fprintf('              Nonstochastic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Asset Dividend   %5.3f         %5.3f         %5.3f\n'  ,[dbar  davg dstd])
fprintf('Asset Price      %5.3f         %5.3f         %5.3f\n\n',[pstar pavg pstd])

% Plot ergodic dividend distribution
[qq,dd] = ksdensity(dsim(:),'support','positive','bandwidth',0.02);
figure
plot(dd,qq)
xlim([0 2])   
xtickformat('%.1f')
title('Ergodic Dividend Distribution')
xlabel('Dividend')
ylabel('Density')

% Plot ergodic asset price distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.02);
figure
plot(pp,qq)
xlim([5 15])  
ytickformat('%.1f')
title('Ergodic Asset Price Distribution')
xlabel('Price')
ylabel('Density')

% Plot ergodic supply distribution
figure
histogram(dsim(:),30,'Normalization','pdf')
xlim([0 2])   
xtickformat('%.1f')
title('Ergodic Dividend Distribution')
xlabel('Dividend')
ylabel('Density')

% Plot ergodic market price distribution
figure
histogram(psim(:),30,'Normalization','pdf')
xlim([5 15])  
ytickformat('%.1f')
title('Ergodic Asset Price Distribution')
xlabel('Price')
ylabel('Density')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
optset('remsolve','nr',0); 
pavgplot = zeros(nplot,1);
pstdplot = zeros(nplot,1);
pinit = p;

fprintf('Varying Coefficient of Risk Aversion\n')
modelplot = model;
alphamin = 0.3;
alphamax = 0.7;
alphaplot = nodeunif(nplot,alphamin,alphamax);
p = pinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {delta alphaplot(ip) dbar gamma};
  [p,c,dr,pr,fr] = remsolve(modelplot,basis,p);
  [dsim,psim] = remsimul(modelplot,basis,nper,dinit,dr,pr,esim);
  pavgplot(ip) = mean(psim(:,nper));
  pstdplot(ip) = std(psim(:,nper));
end
figure
plot(alphaplot,pavgplot)
xticks(alphamin:0.1:alphamax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Coefficient of Risk Aversion')
ylabel('Ergodic Mean Asset Price')
figure
plot(alphaplot,pstdplot)
xticks(alphamin:0.1:alphamax)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Coefficient of Risk Aversion')
ylabel('Ergodic Standard Deviation of Asset Price')


%% SAVE FIGURES
printfigures(mfilename)


%% REMSOLVE FUNCTION FILE
function [out1,out2,out3,out4] = func(flag,d,p,dn,pn,e,delta,alpha,dbar,gamma)
ns = length(d);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1)-inf;
    out2 = zeros(ns,1)+inf;
  case 'f'      % reward
    u = d.^(-alpha);
    un = dn.^(-alpha);
    out1 = p.*u-delta*(pn+dn).*un;
    out2 = u;
    out3 = alpha*delta*(pn+dn).*un./dn-delta*un;
    out4 = -delta*un;
  case 'g'      % transition
    out1 = dbar+gamma*(d-dbar)+e;
    out2 = zeros(ns,1);
end
end