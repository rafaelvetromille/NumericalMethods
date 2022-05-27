%% DEMDP04 Job Search Model
%
% Infinitely-lived worker must decide whether to quit, if employed, or 
% search for a job, if unemployed, given prevailing market wages.

% States
%     w       prevailing wage
%     i       unemployed (1) or employed (2) at beginning of period
% Actions
%     j       idle (1) or active (i.e., work or search) (2) this period
% Parameters
%     l        benefit of pure leisure
%     wbar     long-run mean wage
%     gamma    wage reversion rate
%     p0       probability of finding job
%     p1       probability of keeping job
%     sigma    standard deviation of wage shock
%     delta    discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';                     % basis function type
n     =   151;                          % degree of approximation
wmin  =    60;                          % minimum wage
wmax  =   140;                          % maximum wage
m     =    15;                          % number of wage shocks
nper  =    41;                          % number of periods simulated
nrep  =  1000;                          % number of replications
nplot =   15;                           % number of parameter values plotted

% Base case model parameters
u     =   85;                           % unemployment benefit
l     =   90;                           % benefit of pure leisure
wbar  =  100;                           % long-run mean wage
gamma = 0.40;                           % wage reversion rate
p0    = 0.50;                         	% probability of finding job
p1    = 0.95;                           % probability of keeping job
sigma =    5;                           % standard deviation of wage shock
delta = 0.95;                           % discount factor

% Approximation structure
basis = fundefn(basistype,n,wmin,wmax);

% Continuous wage shock distribution
[e,p] = qnwnorm(m,0,sigma^2);

% Stochastic discrete state transition probabilities
q = zeros(2,2,2);
q(1,2,2) = p0;
q(2,2,2) = p1;
q(:,1,:) = 1-q(:,2,:);

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {u l wbar gamma};      	% function parameters
model.discount = delta;                	% discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 2;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                          	% continuous state shocks
model.w  = p;                       	% continuous state shock probabilities
model.q  = q;                          	% discrete state transition probabilities


%% BASE CASE MODEL SOLUTION

% Solve collocation equation
[v,~,c,wr,vr,~,resid] = dpsolve(model,basis); 

% Critical action wages
wcrit1 = interp1(vr(:,1,1)-vr(:,1,2),wr,0);
vcrit1 = interp1(wr,vr(:,1,1),wcrit1);
wcrit2 = interp1(vr(:,2,1)-vr(:,2,2),wr,0);
vcrit2 = interp1(wr,vr(:,2,1),wcrit2);

% Print critical action wages
fprintf('Critical Wages\n') 
fprintf('   Search  %5.1f\n'  ,wcrit1) 
fprintf('   Quit    %5.1f\n\n',wcrit2)

% Plot action-contingent value functions
figure
hold on
plot(wr,squeeze(vr(:,1,1)))
plot(wr,squeeze(vr(:,1,2)))
plot(wr,squeeze(vr(:,2,2)))
xlim([65 80])
ylim([1920 1945])
plotvdash(wcrit1,vcrit1)
plotbullet(wcrit1,vcrit1)
plotvdash(wcrit2,vcrit2)
plotbullet(wcrit2,vcrit2)
plottext(wcrit1+0.1,[],'$w_0^*$')
plottext(wcrit2+0.1,[],'$w_1^*$')
title('Action-Contingent Value Functions')
xlabel('Wage')
ylabel('Lifetime Utility')
legend('Idle','Search','Work','Location','S')

% Plot residual
figure
hold on
plot(wr,100*resid./max(vr,[],3))
plothdash([],0)
plotvdash(wcrit1,0)
plotbullet(wcrit1,0)
plottext(wcrit1+1,[],'$w_0^*$')
plotvdash(wcrit2,0)
plotbullet(wcrit2,0)
plottext(wcrit2+1,[],'$w_1^*$')
xlim([60 90])
ytickformat('%.0f')
title('Bellman Equation Residual')
xlabel('Wage')
ylabel('Percent Residual')
legend('Unemployed','Employed')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm(0,sigma^2,nrep,nper);

% Initialize simulation
winit = wbar*ones(nrep,1);              % initial wages
iinit = ones(nrep,1);                   % initial employment states

% Simulate model
[wsim,~,isim] = dpsimul(model,basis,nper,winit,iinit,wr,vr,[],esim);

% Ergodic moments
isim = isim-1;
wavg = mean(wsim(:,nper));
iavg = mean(isim(:,nper));
wstd = std(wsim(:,nper));
istd = std(isim(:,nper));
fprintf('                Ergodic       Ergodic\n') 
fprintf('                 Mean      Std Deviation\n') 
fprintf('Wage           %6.2f         %6.2f\n',[wavg wstd])
fprintf('Unemployment   %6.2f         %6.2f\n\n',[1-iavg istd])

% Plot expected discrete state path
figure
hold on
plot(0:nper-1,1-mean(isim),'k')
ytickformat('%.1f')
title('Probability of Unemployment')
xlabel('Period')
ylabel('Probability')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initializations
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
iavgplot = zeros(nplot,1);
vinit = v;

% Sensitivity: Unemployment Benefit
fprintf('Varying Unemployment Benefit\n')
modelplot = model;
umin = 80;
umax = 100;
uplot = nodeunif(nplot,umin,umax);
v = vinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {uplot(ip) l wbar gamma};
  [v,~,c,wr,vr,~] = dpsolve(modelplot,basis);
  [wsim,~,isim] = dpsimul(modelplot,basis,nper,winit,[],wr,vr,[],esim);
  iavgplot(ip) = mean(isim(:,nper));
end
figure
plot(uplot,100*(2-iavgplot))
ytickformat('%.0f%%')
xlabel('Unemployment Benefit')
ylabel('Unemployment Rate')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function out = func(flag,w,~,i,j,e,u,l,wbar,gamma)
switch flag
  case 'f'      % reward
    out = (j==1)*l + (j==2)*((i==1)*u+(i==2)*w);
  case 'g'      % transition
    out = wbar+gamma*(w-wbar)+e;
end
end