%% DEMREM02 Commodity Storage Model
%
% Private expected profit maximizing storers enforce intertemporal price
% equilibrium in a market for a storable commodity.

% States
%     s         supply
% Response
%     x         ending stocks
% Parameters
%     gamma     inverse demand elasticity
%     xbar      storage capacity
%     k         unit storage cost
%     sigma     production volatility
%     delta     discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';                                 % basis function type
n    = 150;                                         % degree of approximation
smin = 0.3;                                         % minimum state
smax = 4.0;                                         % maximum state
m    = 5;                                           % number of shocks
nper = 20;                                          % number of periods simulated
nrep = 10000;                                       % number of replications
nplot = 35;                                         % number of parameter values plotted

% Base case model parameters
gamma = 0.5;                                        % inverse demand elasticity
xbar  = 0.5;                                        % storage capacity
k     = 0.05;                                       % unit storage cost
sigma = 0.3;                                        % production  volatility
delta = 0.95;                                       % discount factor

% Continuous state shock distribution
[y,w] = qnwlogn(m,1,sigma^2);            % shocks and weights

% Model structure
clear model
model.func = @func;                                 % model functions
model.params = {delta,gamma,xbar,k};                % other parameters
model.e = y;                                        % shocks
model.w = w;                                        % shock probabilities

% Approximation structure
[basis,Phi,s] = fundefn(basistype,n,smin,smax);       % basis functions


%% BASE CASE MODEL SOLUTION

% Solve rational expectations equilibrium directly by function iteration
x = zeros(n,1);
c = zeros(n,1);
fprintf('\n')
for it=1:100
  cold = c;
  f = -(s-x).^(-gamma)-k;
  d = -gamma*(s-x).^(-gamma-1);
  for j=1:m
    sn = x + y(j);
    f = f + w(j)*delta*funeval(c,basis,sn);
    d = d + w(j)*delta*funeval(c,basis,sn,1);
  end
  x = x + min(max(-f./d,-x),xbar-x);
  p = (s-x).^(-gamma);
  c = Phi\p;
  change = norm(c-cold,inf);
  fprintf ('%4i %10.1e\n',it,change)
  if change<1.e-10, break, end
end

% Solve rational expectations equilibrium using remsolve
[x,c,sr,xr,fr,resid] = remsolve(model,basis);
pe = (sr-xr).^(-gamma);
p0 = (sr).^(-gamma);

% Critical supply levels
scrit1 = broyden(@(s) delta*w'*((   0+y-funeval(c,basis,   0+y)).^(-gamma))-(s-   0).^(-gamma)-k,1);
scrit2 = broyden(@(s) delta*w'*((xbar+y-funeval(c,basis,xbar+y)).^(-gamma))-(s-xbar).^(-gamma)-k,1);
pcrit1 = interp1(sr,pe,scrit1);
pcrit2 = interp1(sr,pe,scrit2);

% Plot equilibrium ending stocks function
figure
hold on
plot(sr,xr)
plotvdash(scrit1,0)
plotvdash(scrit2,xbar)
plottext(scrit1+0.08,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')
xlim([0 smax])
ylim([0 xbar+0.1])
ytickformat('%.1f')
title('Equilibrium Ending Stocks')
xlabel('Supply')
ylabel('Ending Stocks')

% Plot equilibrium price function
figure
hold on
plot(sr,[pe p0])
plotvdash(scrit1,pcrit1)
plotvdash(scrit2,pcrit2)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')
xlim([0 smax])
ytickformat('%.1f')
title('Equilibrium Market Price')
xlabel('Supply')
ylabel('Price')
legend('Total Demand','Consumption Demand')

% Plot expected arbitrage profit function
figure
hold on
plot(sr,fr)
plothdash([],0)
plotvdash(scrit1,0)
plotvdash(scrit2,0)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')
xlim([0 smax])
ytickformat('%.1f')
title('Expected Arbitrage Profit')
xlabel('Supply')
ylabel('Profit')

% Plot residual
figure
hold on
plot(sr,resid)
plothdash([],0)
plotvdash(scrit1,0)
plotvdash(scrit2,0)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')
xlim([0 smax])
ytickformat('%.1f')
title('Response Residual')
xlabel('Supply')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
ysim = randlogn(1,sigma^2,nrep,nper+1);

% Initialize simulation
sinit = ones(nrep,1);

% Simulate model directly, without remsimul
ssim  = zeros(nrep,nper+1);
xsim  = zeros(nrep,nper+1);
ss = sinit;
for ip=1:nper+1
  xx = interp1(s,x,ss,'PCHIP');
  ssim(:,ip,:) = ss;
  xsim(:,ip,:) = xx;
  if ip<nper+1
    ss = xx + ysim(:,ip);
  end
end

% Simulate model using remsimul
[ssim,xsim] = remsimul(model,basis,nper,sinit,sr,xr,ysim);
psim = (ssim-xsim).^-gamma;

% Plot simulated and expected state path
figure
plot(0:nper,ssim(1:3,:),0:nper,mean(ssim),'k')
ytickformat('%.1f')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot simulated and expected action path
figure
plot(0:nper,xsim(1:3,:),0:nper,mean(xsim),'k')
ytickformat('%.2f')
title('Simulated and Expected Ending Stocks')
xlabel('Period')
ylabel('Ending Stocks')

% Plot simulated and expected market price
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
ytickformat('%.1f')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Market Price')

% Ergodic moments
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 
pstd = std(psim(:)); 

% Print ergodic moments
fprintf('\n')
fprintf('              Nonstochastic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Supply           %5.3f         %5.3f         %5.3f\n',  [1 savg sstd])
fprintf('Ending Stocks    %5.3f         %5.3f         %5.3f\n',  [0 xavg xstd])
fprintf('Market Price     %5.3f         %5.3f         %5.3f\n\n',[1 pavg pstd])

% Plot ergodic supply distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.1);
figure
plot(ss,qq)
xlim([0 3])  
xtickformat('%.1f')
ytickformat('%.1f')
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Density')

% Plot ergodic market price distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.05);
figure
plot(pp,qq)
xlim([0 2])  
xtickformat('%.1f')
ytickformat('%.1f')
title('Ergodic Market Price Distribution')
xlabel('Market Price')
ylabel('Density')

% Plot ergodic supply distribution
figure
histogram(ssim(:),20,'Normalization','pdf')
xlim([0 3])  
xtickformat('%.1f')
ytickformat('%.1f')
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Density')

% Plot ergodic market price distribution
figure
histogram(psim(:),20,'Normalization','pdf')
xlim([0 2])  
xtickformat('%.1f')
ytickformat('%.1f')
title('Ergodic Market Price Distribution')
xlabel('Market Price')
ylabel('Density')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
optset('remsolve','nr',0); 
pavgplot = zeros(nplot,1);
pstdplot = zeros(nplot,1);
xinit = x;

fprintf('Varying Unit Storage Cost\n')
modelplot = model;
kmin = 0.0;
kmax = 0.2;
kplot = nodeunif(nplot,kmin,kmax);
x = xinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {delta gamma xbar kplot(ip)};
  [x,c,sr,xr,fr] = remsolve(modelplot,basis,x);
  [ssim,xsim] = remsimul(modelplot,basis,nper,sinit,sr,xr,ysim);
  psim = (ssim-xsim).^-gamma;
  pavgplot(ip) = mean(psim(:,nper));
  pstdplot(ip) = std(psim(:,nper));
end
figure
plot(kplot,pavgplot)
xticks(kmin:0.1:kmax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Unit Storage Cost')
ylabel('Ergodic Mean Market Price')
figure
plot(kplot,pstdplot)
xticks(kmin:0.1:kmax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Unit Storage Cost')
ylabel('Ergodic Standard Deviation of Market Price')

fprintf('Varying Production Volatility\n')
modelplot = model;
sigmamin = 0.2;
sigmamax = 0.4;
sigmaplot = nodeunif(nplot,sigmamin,sigmamax);
x = xinit;
for ip=1:nplot
  rng('default')
  ysim = randlogn(1,sigmaplot(ip)^2,nrep,nper+1);
  [yplot,wplot] = qnwlogn(m,1,sigmaplot(ip)^2);
  modelplot.e = yplot;
  modelplot.w = wplot;
  [x,c,sr,xr,fr] = remsolve(modelplot,basis,x);
  [ssim,xsim] = remsimul(modelplot,basis,nper,sinit,sr,xr,ysim);
  psim = (ssim-xsim).^-gamma;
  pavgplot(ip) = mean(psim(:,nper));
  pstdplot(ip) = std(psim(:,nper));
end
figure
plot(sigmaplot,pavgplot)
xticks(sigmamin:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Production Volatility')
ylabel('Ergodic Mean Market Price')
figure
plot(sigmaplot,pstdplot)
xticks(sigmamin:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Production Volatility')
ylabel('Ergodic Standard Deviation of Market Price')



%% SAVE FIGURES
printfigures(mfilename)


%% REMSOLVE FUNCTION FILE
function [out1,out2,out3,out4] = func(flag,s,x,sn,xn,y,delta,gamma,xbar,k)
ns = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1);
    out2 = xbar*ones(ns,1);
  case 'f'      % arbitrage profit
    out1 = delta*(sn-xn).^(-gamma)-(s-x).^(-gamma)-k;
    out2 = -gamma*(s-x).^(-gamma-1);
    out3 = -gamma*delta*(sn-xn).^(-gamma-1);
    out4 =  gamma*delta*(sn-xn).^(-gamma-1);
  case 'g'      % transition
    out1 = x + y;
    out2 = ones(ns,1);
end
end