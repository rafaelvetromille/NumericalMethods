%% DEMDP06 Ramsey Stochastic Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest.  

% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     beta	  production elasticity
%     sigma   production shock volatility
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';     % basis function type
n     = 15;             % degree of approximation
smin  = 0.1;            % minimum stock of wealth
smax  = 1.2;            % maximum stock of wealth
m     = 7;              % number of production shocks
nper  = 21;             % number of periods simulated
nrep  = 5000;           % number of replications
nplot = 35;             % number of parameter values plotted
    
% Base case model parameters
beta  = 0.7;            % capital production elasticity
sigma = 0.1;            % production shock volatility
delta = 0.9;            % discount factor

% Approximation structure
[basis,~,snodes] = fundefn(basistype,n,smin,smax);

% Production shock distribution
[e,w] = qnwlogn(m,1,sigma^2);

% Model structure
clear model
model.func = @func;     % model functions
model.params = {beta};  % function parameters
model.discount = delta; % discount factor
model.e = e;            % shocks
model.w = w;            % shock probabilities


%% BASE CASE MODEL SOLUTION

% Nonstochastic steady-state
sstar = (beta*delta)^(beta/(1-beta));   % wealth
kstar = beta*delta*sstar;               % capital investment
vstar = log(sstar-kstar)/(1-delta);     % value
pstar = 1/(sstar*(1-beta*delta));       % shadow price

% Check model derivatives
dpcheck(model,sstar,kstar)

% Closed-form solution
b = 1/(1-delta*beta);
vtrue = @(s) vstar - (0.5*b*delta*sigma^2)/(1-delta) + b*(log(s)-log(sstar));
ktrue = @(s) delta*beta*s;

% Solve collocation equation
[v,k,c,sr,vr,kr,resid] = dpsolve(model,basis,vtrue(snodes),ktrue(snodes));
   
% Linear-Quadratic approximation at refined state grid
[vlq,klq,plq] = lqapprox(model,sr,sstar,kstar);

% Plot optimal policy
figure
hold on
plot(sr,[kr klq])
plothdash(sstar,kstar)
plotvdash(sstar,kstar)
plotbullet(sstar,kstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],kstar,'$k^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Investment')
legend('Chebychev Collocation','L-Q Approximation')

% Plot value function
figure
hold on
plot(sr,[vr vlq])
plotvdash(sstar,vstar)
plothdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],vstar,'$v^*$')
xtickformat('%.1f')
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')
legend('Chebychev Collocation','L-Q Approximation')

% Plot shadow price function
figure
hold on
pr = funeval(c,basis,sr,1);
plot(sr,[pr plq])
plothdash(sstar,pstar)
plotvdash(sstar,pstar)
plotbullet(sstar,pstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],pstar-0.5,'$\lambda^*$')
xtickformat('%.1f')
title('Shadow Price Function')
xlabel('Wealth')
ylabel('Shadow Price')
legend('Chebychev Collocation','L-Q Approximation')

% Plot collocation equation residual and Bellman equation approximation error
figure
hold on
plot(sr,[resid vr-vtrue(sr)])
plothdash([],0)
xtickformat('%.1f')
ytickformat('%.1f')
title('Collocation Equation Residual and Belman Equation Approximation Error')
xlabel('Wealth')
ylabel('Residual/Error')
legend('Residual','Error')

% Plot linear-quadratic Belman equation approximation error
figure
hold on
plot(sr,vlq-vtrue(sr))
xtickformat('%.1f')
ytickformat('%.1f')
title('Linear-Quadratic Bellman Equation Approximation Error')
xlabel('Wealth')
ylabel('Error')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(1,sigma^2,nrep,nper);

% Initialize simulation
sinit = smin*ones(nrep,1);

% Simulate model
[ssim,ksim] = dpsimul(model,basis,nper,sinit,[],sr,vr,kr,esim);

% Plot simulated and expected wealth path
figure
hold on
if sigma>0
  plot(0:nper-1,ssim(1:3,:))
  plot(0:nper-1,mean(ssim),'k')
else
  plot(0:nper-1,ssim)
end
plothdash([],sstar)
plottext([],sstar,'$s^*$')
ytickformat('%.1f')
title('Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

% Plot simulated and expected investment path
figure
hold on
if sigma>0
  plot(0:nper-1,ksim(1:3,:))
  plot(0:nper-1,mean(ksim),'k')
else
  plot(0:nper-1,ksim)
end
plothdash([],kstar)
plottext([],kstar,'$k^*$')
ytickformat('%.2f')
title('Simulated and Expected Investment')
xlabel('Period')
ylabel('Investment')

% Ergodic moments
ssim = ssim(:,nper);
ksim = ksim(:,nper);
savg = mean(ssim(:));
kavg = mean(ksim(:));
sstd = std(ssim(:));
kstd = std(ksim(:));

% Print ergodic moments
fprintf('          Nonstochasic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Wealth       %5.3f         %5.3f         %5.3f\n'  ,[sstar savg sstd])
fprintf('Investment   %5.3f         %5.3f         %5.3f\n\n',[kstar kavg kstd])

% Plot ergodic wealth distribution
if sigma>0
  [qq,ss] = ksdensity(ssim(:),'support','positive');
  figure
  plot(ss,qq)
  xlim([0.1 0.6])
  title('Ergodic Distribution of Wealth')
  xlabel('Wealth')
  ylabel('Probability Density')
end

% Plot ergodic wealth distribution
if sigma>0
  figure
  histogram(ssim(:),'Normalization','pdf')
  xlim([0.1 0.6])
  title('Ergodic Distribution of Wealth')
  xlabel('Wealth')
  ylabel('Probability Density')
end


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
savgplot = zeros(nplot,1);
sstdplot = zeros(nplot,1);
vinit = v;
kinit = k;

fprintf('Varying Production Elasticity\n')
modelplot = model;
betamin = 0.6;
betamax = 0.8;
betaplot = nodeunif(nplot,betamin,betamax);
for ip=1:nplot
  modelplot.params = {betaplot(ip)};
  sstar = (betaplot(ip)*delta)^(betaplot(ip)/(1-betaplot(ip)));
  kstar = betaplot(ip)*delta*sstar;
  vstar = log(sstar-kstar)/(1-delta);
  b = 1/(1-delta*betaplot(ip));
  v = vstar - (0.5*b*delta*sigma^2)/(1-delta) + b*(log(snodes)-log(sstar));
  k = delta*betaplot(ip)*snodes;
  [~,~,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(betaplot,savgplot)
xlabel('Production Elasticity')
ylabel('Ergodic Mean Wealth')
xlim([betamin betamax])
xticks(betamin:0.1:betamax)
ytickformat('%.2f')
figure
plot(betaplot,sstdplot)
xlabel('Production Elasticity')
ylabel('Ergodic Standard Deviation of Wealth')
xlim([betamin betamax])
xticks(betamin:0.1:betamax)
ytickformat('%.3f')

fprintf('Varying Discount Factor\n')
modelplot = model;
deltamin = 0.85;
deltamax = 0.95;
deltaplot = nodeunif(nplot,deltamin,deltamax);
for ip=1:nplot
  modelplot.discount = deltaplot(ip);
  sstar = (beta*deltaplot(ip))^(beta/(1-beta));
  kstar = beta*deltaplot(ip)*sstar;
  vstar = log(sstar-kstar)/(1-deltaplot(ip));
  b = 1/(1-deltaplot(ip)*beta);
  v = vstar - (0.5*b*deltaplot(ip)*sigma^2)/(1-deltaplot(ip)) + b*(log(snodes)-log(sstar));
  k = deltaplot(ip)*beta*snodes;
  [~,~,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(deltaplot,savgplot)
xlabel('Discount Factor')
ylabel('Ergodic Mean Wealth')
xtickformat('%.2f')
ytickformat('%.2f')
figure
plot(deltaplot,sstdplot)
xlabel('Discount Factor')
ylabel('Ergodic Standard Deviation of Wealth')
xtickformat('%.2f')
ytickformat('%.3f')

fprintf('Varying Production Shock Volatility\n')
modelplot = model;
sigmamin = 0.0;
sigmamax = 0.2;
sigmaplot = nodeunif(nplot,sigmamin,sigmamax);
v = vinit; k = kinit;
for ip=1:nplot
  [e,w] = qnwlogn(m,1,sigmaplot(ip)^2);
  modelplot.e = e;
  modelplot.w = w;
  rng('default')
  esim = randlogn(1,sigmaplot(ip)^2,nrep,nper);
  [v,k,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(sigmaplot,savgplot)
xticks(0:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Production Shock Volatility')
ylabel('Ergodic Mean Wealth')
figure
plot(sigmaplot,sstdplot)
xticks(0:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Production Shock Volatility')
ylabel('Ergodic Standard Deviation of Wealth')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,k,~,~,e,beta)
n = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(n,1);
    out2 = s;
    out3 = [];
  case 'f'      % reward  
    out1 = log(s-k);
    out2 = -(s-k).^(-1);
    out3 = -(s-k).^(-2);
  case 'g'      % transition
    out1 = e.*k.^beta;
    out2 = beta*e.*k.^(beta-1);
    out3 = (beta-1)*beta*e.*k.^(beta-2);
end
end