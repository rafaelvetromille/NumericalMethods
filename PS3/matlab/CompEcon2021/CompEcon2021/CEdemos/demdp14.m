%% DEMDP14 Livestock Feeding Model
%
% Farmer must decide how much to feed his livestock over a finite horizon.

%  States
%      s       weight of livestock at beginning of period
%  Actions
%      x       weight gain current period
%  Parameters
%      alpha   weight gain function parameter
%      beta    weight gain function parameter
%      k       unit cost of feed
%      p       price of livestock per pound
%      N       number of feeding periods
%      s1      initial livestock weight
%      delta   discount factor

% Preliminary tasks
deminit(mfilename)



%% FORMULATION

% Numerical control parameters
basistype = 'spli';                    % basis function type
n     = 50;                            % degree of approximation
smax  = 4.0;                           % maximum livestock weight
nplot =  15;                           % number of parameter values plotted
 
% Base case model parameters
alpha = 0.1;                           % weight gain function parameter
beta  = 2.0;                           % weight gain function parameter
k     = 0.4;                           % unit cost of feed
p     = 1.0;                           % price of livestock per pound
T     = 6;                             % number of feeding periods
s1    = 0.4;                           % initial livestock weight
delta = 0.95;                          % per-period discount factor

% Approximation structure
[basis,~,snodes] = fundefn(basistype,n,s1,smax);

% Model structure
clear model
model.horizon = T;                     % number of decision periods
model.func = @func;                    % model function
model.discount = delta;                % discount factor
model.params = {alpha beta k};         % other parameters

% Post-terminal value function
vterm = p*snodes;


%% BASE CASE MODEL SOLUTION

% Solve collocation equation
[v,x,c,sr,vr,xr] = dpsolve(model,basis,vterm);

% Plot optimal policy
figure
plot(sr,xr)
ytickformat('%.1f')
title('Optimal Weight Gain Policy')
xlabel('Livestock Weight')
ylabel('Weight Gain')
legend('t=1','t=2','t=3','t=4','t=5','t=6')

% Plot value function
figure
plot(sr,vr)
title('Value Function')
xlabel('Livestock Weight')
ylabel('Present Value of Current and Future Profit')
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7')

% Plot shadow price function
figure
p = funeval(c,basis,sr,1);
plot(sr,p)
ytickformat('%.1f')
title('Shadow Price Function')
xlabel('Livestock Weight')
ylabel('Shadow Price')
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7')


%% BASE CASE MODEL SIMULATION

% Simulate model
[ssim,xsim] = dpsimul(model,basis,T,s1,1,sr,vr,xr);

% Terminal weight
fprintf('\n')
fprintf('Terminal Weight = %5.2f\n',ssim(end))

% Plot simulated state path
figure
plot(1:T+1,ssim)
xlim([1 T+1])
title('Simulated Livestock Weight')
xlabel('Period')
ylabel('Weight')

% Plot simulated action path
figure
plot(1:T,xsim)
xlim([1 T])
ytickformat('%.2f')
title('Simulated Livestock Weight Gain')
xlabel('Period')
ylabel('Weight Gain')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
termweightplot = zeros(nplot,1);

fprintf('Varying Unit Cost of Feed\n')
modelplot = model;
kmin = 0.1;
kmax = 0.6;
kplot = nodeunif(nplot,kmin,kmax);
for ip=1:nplot
  modelplot.params = {alpha beta kplot(ip)};
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,vterm);
  ssim = dpsimul(modelplot,basis,T,s1,1,sr,vr,xr);  
  termweightplot(ip) = ssim(end);
end
figure
plot(kplot,termweightplot)
xlabel('Unit Cost of Feed')
ylabel('Terminal Weight')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,x,~,~,~,alpha,beta,k)
ns = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1);
    out2 = inf*ones(ns,1);
    out3 = [];
  case 'f'      % reward
    out1 = -k*(x+alpha*s).^beta;
    out2 = -k*beta*(x+alpha*s).^(beta-1);
    out3 = -k*beta*(beta-1)*(x+alpha*s).^(beta-2);
  case 'g'      % transition
    out1 = s + x;
    out2 = ones(ns,1);
    out3 = zeros(ns,1);
end
end