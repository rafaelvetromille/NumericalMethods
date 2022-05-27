%% DEMSOC03 Stochastic Optimal Economic Growth Model
%
% Social benefit maximizing social planner must decide how much society
% should consume and invest.

% State
%     k       capital stock
%     y       productivity shock
% Control
%     q       consumption rate
% Parameters
%     alpha   capital share
%     delta   capital depreciation rate
%     theta   relative risk aversion
%     gamma   productivity mean reversion coefficient
%     sigma   productivity volatility
%     rho     discount rate

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';     % basis function type
n     = [15 5];         % degree of approximation
kmin  =  2;             % minimum capital stock
kmax  =  8;             % maximum capital stock
ymin  = 0.9;            % minimum productivity factor
ymax  = 1.1;            % maximum productivity factor
sinit = [3 1];          % initial capital stock and productivity shock
T     = 25;             % time horizon
nrep  = 2000;           % number of replications

% Base case model parameters
alpha =  0.4;                                   % capital share
delta =  0.1;                                   % capital depreciation rate
theta =  2.0;                                   % relative risk aversion
gamma =  0.5;                                   % productivity mean reversion coefficient
sigma =  0.05;                                  % productivity volatility
rho   =  0.05;                                  % discount rate

% Model Structure
clear model
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha delta gamma theta sigma}; % function file parameters

% Approximation Structure
smin = [kmin ymin];                             % minimum states
smax = [kmax ymax];                             % maximum states
[basis,Phi,snode] = fundefn(basistype,n,smin,smax);  % basis functions


%% BASE CASE MODEL SOLUTION

% Solve HJB equation by collocation
v = (((rho*snode(:,1)).^(1-theta))/(1-theta));
[c,s,v,q,res] = socsolve(model,basis,v);

% Fix y=1
p = funeval(c,basis,s,[1 0]);
k = s(:,1);
y = s(:,2);
j = find(y==1);
k = k(j);
v = v(j);
q = q(j);
p = p(j);
res = res(j);

% Plot optimal policy
figure
plot(k,q)
ytickformat('%.1f')
title('Optimal Consumption Policy')
xlabel('Capital Stock')
ylabel('Rate of Consumption')

% Plot value function
figure
plot(k,v)
title('Value Function')
xlabel('Capital Stock')
ylabel('Lifetime Utility')

% Plot shadow price function
figure
plot(k,p)
ytickformat('%.1f')
title('Shadow Price Function')
xlabel('Capital Stock')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(k,res)
plothdash([],0)
title('HJB Equation Residual')
xlabel('Capital Stock')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Simulate model
[t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);
ksim = ssim(:,:,1);
ysim = ssim(:,:,2);
kmean = smean(:,1);
ymean = smean(:,2);

% Plot simulated and expected state path
figure
plot(t,ksim,t,kmean,'k')
ytickformat('%.1f')
title('Simulated and Expected Capital Stock')
xlabel('Time')
ylabel('Capital Stock')

% Plot simulated and expected state path
figure
plot(t,ysim,t,ymean,'k')
ytickformat('%.2f')
title('Simulated and Expected Productivity Shock')
xlabel('Time')
ylabel('Productivity Shock')

% Plot expected control path
figure
plot(t,xsim,t,xmean,'k')
ytickformat('%.1f')
title('Simulated and Expected Rate of Consumption')
xlabel('Time')
ylabel('Rate of Consumption')


%% SOLVE HJB EQUATION DIRECTLY USING BROYDEN

% Define Residual Function
k     = @(s) s(:,1);
y     = @(s) s(:,2);
V     = @(c,s) funeval(c,basis,s);
Vk    = @(c,s) funeval(c,basis,s,[1 0]);
Vy    = @(c,s) funeval(c,basis,s,[0 1]);
Vyy   = @(c,s) funeval(c,basis,s,[2 2]);
q     = @(c,s) Vk(c,s).^(-1/theta);
resid = @(c,s) (q(c,s).^(1-theta))/(1-theta) ...
       + (y(s).*k(s).^alpha-delta*k(s)-q(c,s)).*Vk(c,s) ...
       + gamma*(1-y(s)).*Vy(c,s) ...
       - 0.5*(sigma^2)*y(s).*Vyy(c,s) - rho*V(c,s);

% Solve
v = (((rho*k(snode)).^(1-theta))/(1-theta));
c = Phi\v;
c = broyden(resid,c,snode);
% norm(resid(c,s))


%% SAVE FIGURES
printfigures(mfilename)


%% SOCSOLVE FUNCTION FILE
function out = func(flag,s,q,Vs,~,alpha,delta,gamma,theta,sigma)
k = s(:,1);
y = s(:,2);
switch flag
  case 'x'        % optimal control
    Vk  = Vs(:,1);
    out = Vk.^(-1/theta);
  case 'f'        % reward
    out = (1/(1-theta))*q.^(1-theta);
  case 'mu'       % state drift
    f = k.^alpha;
    out = [(y.*f-delta*k-q)  gamma*(1-y)];
  case 'sigma'    % state diffusion
    n = length(k);
    out = zeros(n,2,2);
    out(:,2,2) = sigma*sqrt(y);
end
end