%% DEMDP03 Dixit Industry Entry-Exit Model
%
% Profit maximizing firm must decide whether to operate or shut down, given
% its current operational state and current profitability, subject to
% transactions costs.

% States
%     p       current profitability
%     i       active (1) or idle (2) last period
% Actions
%     j       active (1) or idle (2) this period
% Parameters
%     pbar    long-run mean profitability
%     gamma   profitability autoregressive coefficient
%     kappa   cost of reopenning idle firm
%     sigma   standard deviation of profitability shock
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb'; % basis function type
n     =  150;       % degree of approximation
pmin  =  -15;       % minimum price state
pmax  =   15;       % maximum price state
m     =   21;       % number of price shocks
nper  =   21;       % number of periods simulated
nrep  = 10000;      % number of replications simulated
nplot =   15;       % number of parameter values plotted
  
% Base case model parameters
pbar  = 1.0;        % long-run mean profitability
gamma = 0.7;        % profitability autocorrelation
kappa = 1.0;        % cost of reopenning idle firm
sigma = 1.0;        % standard deviation of profitability shock
delta = 0.9;        % discount factor

% Approximation structure
basis = fundefn(basistype,n,pmin,pmax);
  
% Profitability innovation distribution
[e,w] = qnwnormunif(m,0,sigma.^2);
  
% Model structure
clear model
model.func = @func;                     % model functions
model.params = {pbar gamma kappa};      % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 2;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities
model.h  = [1 1;2 2];                   % deterministic discrete state transitions


%% BASE CASE MODEL SOLUTION

% Solve collocation equation
[v,~,~,pr,vr,~,resid] = dpsolve(model,basis); 

% Critical profitabilities and values
pcrit = [interp1(vr(:,1,1)-vr(:,1,2),pr,0); interp1(vr(:,2,1)-vr(:,2,2),pr,0)];
vcrit = [interp1(pr,vr(:,1,1),pcrit(1)); interp1(pr,vr(:,2,1),pcrit(2))];

% Print critical profitabilities
fprintf('Critical Profitabilities\n')
fprintf(' Exit  %5.2f\n'  ,pcrit(1))  
fprintf(' Entry %5.2f\n\n',pcrit(2))

% Plot action-contingent value functions
figure
hold on
plot(pr,[vr(:,1,1) vr(:,2,1) vr(:,2,2)])
for i=1:2
  plotvdash(pcrit(i),vcrit(i))
  plotbullet(pcrit(i),vcrit(i))
end
plottext(pcrit(1)+0.05,0,'$p^*_1$')
plottext(pcrit(2)+0.05,0,'$p^*_0$')
xlim([-1 1])
ylim([0 15])
xtickformat('%.1f')
title('Action-Contingent Value Functions')
xlabel('Profitability')
ylabel('Value of Firm')
legend('Keep Active Firm Open','Reopen Idle Firm','Shut Down')

% Plot residual
figure
hold on
plot(pr,100*resid./max(vr,[],3))
plothdash([],0)
for i=1:2
  plotvdash(pcrit(i),0)
  plotbullet(pcrit(i),0)
end
plottext(pcrit(1)+0.1,[],'$p^*_1$')
plottext(pcrit(2)+0.1,[],'$p^*_0$')
xlim([-2 2])
ytickformat('%.1f%%')
title('Bellman Equation Residual')
xlabel('Profitability')
ylabel('Percent Residual')
legend('Active','Idle')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm(0,sigma^2,nrep,nper);

% Initialize simulation
pinit = pbar*ones(nrep,1);              % initial profitabilities
iinit = ones(nrep,1);                   % initial firm activity states

% Simulate model
[psim,~,isim,jsim] = dpsimul(model,basis,nper,pinit,iinit,pr,vr,[],esim);

% Convert to 0=idle, 1=active
isim = 2 - isim; 
jsim = 2 - jsim; 

% Compute profit
psim = jsim.*psim - kappa.*jsim.*(1-isim);

% Plot simulated and expected profitability path
figure
hold on
plot(0:nper-1,psim(1:3,:))
plot(0:nper-1,mean(psim),'k')
title('Simulated and Expected Profitabilities')
xlabel('Period')
ylabel('Profitability')

% Plot expected activity path
figure
hold on
plot(0:nper-1,mean(jsim),'k')
ytickformat('%.2f')
title('Probability of Operation')
xlabel('Period')
ylabel('Probability')

% Plot simulated and expected net profit
figure
hold on
plot(0:nper-1,zeros(nper,1),':k')
plot(0:nper-1,psim(1:3,:))
plot(0:nper-1,mean(psim),'k')
title('Simulated and Expected Net Profit')
xlabel('Period')
ylabel('Net Profit')

% Ergodic moments
pavg = mean(psim(:));
pstd = std(psim(:));
javg = mean(jsim(:));
jstd = std(jsim(:));

% Print ergodic moments
fprintf('                Ergodic       Ergodic\n') 
fprintf('                 Mean      Std Deviation\n') 
fprintf('Profit          %5.3f         %5.3f\n'  ,[pavg pstd])
fprintf('Activity        %5.3f         %5.3f\n\n',[javg jstd])


%% PARAMETRIC SENSITIVITY ANALYSIS

% Intialization
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
savgplot = zeros(nplot,1);
iavgplot = zeros(nplot,1);

fprintf('Varying Cost of Reopenning\n')
modelplot = model;
kappamin = 0.0;
kappamax = 3.0;
kappaplot = nodeunif(nplot,kappamin,kappamax);
for ip=1:nplot
  rng('default')
  modelplot.params = {pbar gamma kappaplot(ip)}; 
  [~,~,~,pr,vr,~] = dpsolve(modelplot,basis);
  [psim,~,isim] = dpsimul(modelplot,basis,nper,pinit,iinit,pr,vr,[],esim);
  savgplot(ip) = mean(psim(:,nper));
  iavgplot(ip) = 100*(2-mean(isim(:,nper)));
end
figure
plot(kappaplot,iavgplot)
xticks(0:1:3)
xtickformat('%.0f')
ytickformat('%.0f%%')
xlabel('Cost of Reopenning')
ylabel('Ergodic Percent of Time Active')

fprintf('Varying Profitability Autocorrelation\n')
modelplot = model;
gammamin = 0.0;
gammamax = 0.8;
gammaplot = nodeunif(nplot,gammamin,gammamax);
for ip=1:nplot
  modelplot.params = {pbar gammaplot(ip) kappa}; 
  [~,~,~,pr,vr,~] = dpsolve(modelplot,basis);
  [psim,~,isim] = dpsimul(modelplot,basis,nper,pinit,iinit,pr,vr,[],esim);
  savgplot(ip) = mean(psim(:,nper));
  iavgplot(ip) = 100*(2-mean(isim(:,nper)));
end
figure
plot(gammaplot,iavgplot)
xtickformat('%.1f')
ytickformat('%.0f%%')
xlabel('Profitability Autocorrelation')
ylabel('Ergodic Percent of Time Active')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function out = func(flag,p,~,i,j,e,pbar,gamma,kappa)
switch flag
  case 'f'      % reward
    out = p.*(j==1) - kappa.*(i==2).*(j==1);
  case 'g'      % transition
    out = pbar+gamma*(p-pbar)+e;
end
end