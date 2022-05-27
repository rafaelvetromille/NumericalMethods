%% DEMDP11 Monetary Policy Model
%
% A central bank must set nominal interest rate so as to minimize
% deviations of inflation rate and GDP gap from established targets.

% States
%     s1      GDP gap
%     s2      inflation rate
% Actions
%     x       nominal interest rate
% Parameters
%     alpha   transition function constant coefficients
%     beta    transition function state coefficients
%     gamma   transition function action coefficients
%     omega   central banker's preference weights
%     sbar    equilibrium targets
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)



%% FORMULATION

% Numerical control parameters
basistype = 'spli';                     % basis function type
n = [30 30];                            % degree of approximation
smin = [-2 -3];                         % minimum states
smax = [ 2  3];                         % maximum states
nper = 31;                              % number of periods simulated
nrep = 2000;                            % number of replications
 
% Base case model parameters
alpha   = [0.9; -0.1];                  % transition function constant coefficients
beta    = [-0.5 0.2;0.3 -0.4];          % transition function state coefficients
gamma   = [-0.1; 0.0];                  % transition function action coefficients
omega   = [1 0;0 1];                    % central banker's preference weights
sbar    = [1;0];                        % equilibrium targets
sigma   = 0.2*eye(2);                   % shock covariance matrix
delta   = 0.9;                          % discount factor

% Approximation structure
basis = fundefn(basistype,n,smin,smax);
 
% Continuous state shock distribution
m   = [3 3];                            % number of shocks
[e,w] = qnwnorm(m,[0 0],sigma);         % shocks and probabilities

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {alpha beta gamma omega sbar};	% function parameters
model.discount = delta;                 % discount factor
model.ds = 2;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities


%% BASE CASE MODEL SOLUTION

% Unconstrained nonstochastic steady-state
F0  = -0.5*sbar'*omega*sbar;
Fs  = sbar'*omega;
Fx  = 0;
Fss = -omega;
Fsx = [0 0]';
Fxx = 0;
G0  = alpha;
Gs  = beta;
Gx  = gamma;
[~,~,~,sstar,xstar] = lqsolve(F0,Fs,Fx,Fss,Fsx,Fxx,G0,Gs,Gx,delta);

% Re-compute nonstochastic steady-state, if nonnegativity violated 
if xstar<0
  I = eye(2,2);
  xstar = 0;
  sstar = (I-beta)\alpha;
end

% Reorient nonstochastic steady state
sstar = sstar';

% Check model derivatives
dpcheck(model,sstar,xstar)

% Solve collocation equation
optset('dpsolve','nr',5)
[v,x,c,sr,vr,xr,resid] = dpsolve(model,basis);
optset('dpsolve','nr',10)

% Reshape output for plotting
n = n*5+1;
s1 = reshape(sr(:,1),n);
s2 = reshape(sr(:,2),n);
vr = reshape(vr,n);
xr = reshape(xr,n);
resid = reshape(resid,n);

% Plot optimal policy
figure
mesh(s1,s2,xr)
ytickformat('percentage')
ztickformat('percentage')
title('Optimal Monetary Policy')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Nominal Interest Rate')
  
% Plot value function
figure
mesh(s1,s2,vr)
ytickformat('percentage')
title('Value Function')
xlabel('GDP Gap');
ylabel('Inflation Rate')
zlabel('Social Welfare')

% Plot shadow price function 1
p1 = funeval(c,basis,sr,[1 0]);
p1 = reshape(p1,n);
figure
mesh(s1,s2,p1)
ytickformat('percentage')
title('Shadow Price of GDP Gap')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Price')

% Plot shadow price function 2
p2 = funeval(c,basis,sr,[0 1]);
p2 = reshape(p2,n);
figure
mesh(s1,s2,p2)
ytickformat('percentage')
title('Shadow Price of Inflation Rate')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Price')

% Plot residual
figure
mesh(s1,s2,resid)
ytickformat('percentage')
title('Bellman Equation Residual')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm([0 0],sigma,nrep,nper);

% Initialize simulation
sinit = smax(ones(nrep,1),:);           % initial GDP gap and inflation rate

% Simulate model
[ssim,xsim] = dpsimul(model,basis,nper,sinit,[],sr,vr,xr,esim);
s1sim = ssim(:,:,1);
s2sim = ssim(:,:,2);

% Ergodic moments
s1avg = mean(s1sim(:)); 
s2avg = mean(s2sim(:)); 
xavg  = mean(xsim(:)); 
s1std = std(s1sim(:)); 
s2std = std(s2sim(:)); 
xstd  = std(xsim(:)); 

% Print ergodic moments
fprintf('\n') 
fprintf('Ergodic Moments\n') 
fprintf('                     Nonstochastic    Ergodic      Ergodic\n') 
fprintf('                     Steady-State      Mean     Std Deviation\n') 
fprintf('GDP Gap                 %5.3f         %5.3f         %5.3f\n',[sstar(1) s1avg s1std])
fprintf('Inflation Rate          %5.3f         %5.3f         %5.3f\n',[sstar(2) s2avg s2std])
fprintf('Nominal Interest Rate   %5.3f         %5.3f         %5.3f\n',[xstar    xavg  xstd])

% Plot simulated state path 1
figure
hold on
plot(0:nper-1,s1sim(1:3,:))
plot(0:nper-1,mean(s1sim),'k')
ytickformat('%.1f')
title('Simulated and Expected GDP Gap')
xlabel('Period')
ylabel('GDP Gap')

% Plot simulated state path 2
figure
hold on
plot(0:nper-1,s2sim(1:3,:))
plot(0:nper-1,mean(s2sim),'k')
ytickformat('percentage')
title('Simulated and Expected Inflation Rate')
xlabel('Period')
ylabel('Inflation Rate')

% Plot simulated action path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
plot(nper-1,xavg,'k*')
ytickformat('%.1f%%')
title('Simulated and Expected Nominal Interest Rate')
xlabel('Period')
ylabel('Nominal Interest Rate')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,x,~,~,e,alpha,beta,gamma,omega,sbar)
[n, ds]  = size(s);
dx = 1;
switch flag
  case 'b'      % bounds
    out1  = -0*ones(n,1);
    out2  = inf*ones(n,1);
    out3 = [];
  case 'f'      % reward
    s = s - sbar(:,ones(1,n))';
    out1 = zeros(n,1);
    for i=1:2
      for j=1:2
        out1 = out1 - 0.5*omega(i,j)*s(:,i).*s(:,j);
      end
    end
    out2 = zeros(n,1);
    out3 = zeros(n,1);
  case 'g'      % transition
    out1 = alpha(:,ones(1,n))' + s*beta' + x*gamma' + e;
    out2 = gamma(:,ones(1,n))';
    out3 = zeros(n,ds,dx,dx);
    out2 = reshape(out2,n,ds,dx);
end
end