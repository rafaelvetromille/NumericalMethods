%% DEMDP09 Private Non-Renewable Resource Model
%
% Profit maximizing mine owner must decide how much ore to extract.

% States
%     s       ore stock
% Actions
%     q       ore extracted and sold
% Parameters
%     a       demand function parameters
%     b       cost function parameters
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';                     % basis function type
n    = 101;                             % degree of approximation
smin =   0;                             % minimum ore stock
smax =  10;                             % maximum ore stock
nper =  21;                             % number of periods simulated
nrep =   1;                             % number of replications
 
% Base case model parameters
a = [5 0.8];                            % demand function parameters
b = [7 1.0];                            % cost function parameters
delta = 0.9;                            % discount factor

% Approximation structure
[basis,~,snodes] = fundefn(basistype,n,smin,smax);

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {a b};                   % function parameters
model.discount = delta;                 % discount factor


%% BASE CASE MODEL SOLUTION

% Check model derivatives
dpcheck(model,smax,0)

% Solve collocation equation
[v,q,c,sr,vr,qr,resid] = dpsolve(model,basis);

% Abandonment point
sstar = (b(1)-a(1))/b(2);
fprintf('\n')
fprintf('Abandonment Point %5.2f\n',sstar)

% Plot optimal policy
figure
hold on
plot(sr,qr)
plotvdash(sstar,0)
plotbullet(sstar,0)
plottext(sstar-0.1,[],'$s^*$')
ytickformat('%.1f')
title('Optimal Extraction Policy')
xlabel('Ore Stock')
ylabel('Quantity Extracted')

% Plot value function
figure
plot(sr,vr)
plotbullet(sstar,0)
plottext(sstar-0.1,0,'$s^*$')
ylim([0 20])
title('Value Function')
xlabel('Ore Stock')
ylabel('Value of Mine')

% Plot shadow price function
figure
hold on
pr = funeval(c,basis,sr,1);
plot(sr,pr)
plotvdash(sstar,0)
plotbullet(sstar,0)
plottext(sstar-0.1,0,'$s^*$')
ylim([0 5])
title('Shadow Price Function')
xlabel('Ore Stock')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(sr,resid)
plotbullet(sstar,0)
plottext(sstar+0.1,0,'$s^*$')
ytickformat('%.1f')
title('Bellman Equation Residual')
xlabel('Ore Stock')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Initialize simulation
sinit = smax;

% Simulate model
[ssim,qsim] = dpsimul(model,basis,nper,sinit,[],sr,vr,qr);

% Plot simulated state and policy path
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
plothdash([],sstar)
plottext([],sstar+0.01,'$s^*$')
title('Simulated Ore Stock and Quantity Extracted')
xlabel('Period')
ylabel('Quantity')
legend('Ore Stock','Quantity Extracted')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,q,~,~,~,a,b)
switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (a(1)-b(1)+b(2)*s).*q - (a(2)+b(2)/2).*q.^2;
    out2 = (a(1)-b(1)+b(2)*s) - 2*(a(2)+b(2)/2).*q;
    out3 = -2*(a(2)+b(2)/2)*ones(size(s));
  case 'g'      % transition
    out1 = s-q;
    out2 = -ones(size(s));
    out3 = zeros(size(s));
end
end