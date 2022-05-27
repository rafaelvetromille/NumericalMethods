%% DEMDP20 Saving with Transactions Costs Model
%
% An infinitely-lived agent subject to income uncertainty must decide how
% much to consume, borrow, and save, given that he incurs a fixed
% transaction cost whenever he adjusts his assets and given interest rates
% on debt and deposits differ.

% States
%      s        beginning assets (s<0 indicates indebtedness)
%      i        income state
% Actions
%      send     ending assets
%      j        transaction decision (no=1, yes=2)

% Preliminary tasks
deminit(mfilename)
optset('broyden','maxit',500)

   
%% FORMULATION

% Numerical control parameters
nper = 500;             % number of periods simulated
nrep = 4000;            % number of replications
nwrm = 100;             % number of `warm-up' periods

% Base case model parameters
delta = 0.95;           % discount factor
alpha = 2.0;            % relative risk aversion
smin  = -0.5;           % minimum assets (debt limit)
smax  =  1.0;           % maximum assets (deposit limit)
y     = [1.0;0.5];      % high-low incomes
w     = [0.7;0.3];      % probabilities of high-low incomes
ny    = length(y);      % number of income states

for icase=1:3
  
if icase==1
  tau   = 0.00;         % fixed transactions cost
  rb    = 0.05;         % interest rate on debt
  rs    = 0.05;         % interest rate on deposits
end
if icase==2
  tau   = 0.02;         % fixed transactions cost
  rb    = 0.05;         % interest rate on debt
  rs    = 0.05;         % interest rate on deposits
end
if icase==3
  tau   = 0.00;         % fixed transactions cost
  rb    = 0.05;         % interest rate on debt
  rs    = 0.15;         % interest rate on deposits
end

% Check borrowing limit
if min(y)+rb*smin<=0, warning('Inadmissible borrowing limit.'), end

% Utility function
util = @(c) c.^(1-alpha)/(1-alpha);

% Discretize admissible asset holdings
n = 601;
s = nodeunif(n,smin,smax);
s = unique([s;0]);
n = length(s);

% State matrix for discrete optimization
S = s(:,ones(1,n));

% Net interest
Is = rb*min(s,0) + rs*max(s,0);
IS = rb*min(S,0) + rs*max(S,0);


%% BASE CASE MODEL SOLUTION

% Intialize
v  = zeros(n,ny);
v1 = zeros(n,ny);
v2 = zeros(n,ny);

% Perform some function iterations on collocation equations
fprintf('\n')
for it=1:500
  vold = v;
  Evnext1 = v*w;
  Evnext2 = Evnext1(:,ones(1,n))';
  for i=1:ny
    % do not transact
    v1(:,i) = util(y(i)+Is) + delta*Evnext1;
    % transact
    C = y(i) + IS + S - S' - tau;
    U = util(C);
    U(C<=0) = -inf;
    v2(:,i) = max(U+delta*Evnext2,[],2);
  end
  v = max(v1,v2);
  change = norm(v(:)-vold(:));
%   fprintf ('%4i %10.1e\n',it,change)
  if change<1.e-8, break, end
end
if it==500, warning('demdp15: Failure to Converge.'), end

% Recover optimal policy
[r,send] = resid(v,n,ny,util,s,S,Is,IS,y,w,tau,delta);

% Plot optimal policy
figure
hold on
plot(s,send)
plot(s,s,'k:','LineWidth',2)
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Saving Policy')
xlabel('Beginning Assets')
ylabel('Ending Assets')
legend('High Income','Low Income') 

% Plot value function
figure
hold on
v = reshape(v,n,ny);
plot(s,v)
xtickformat('%.1f')
title('Value Function')
xlabel('Beginning Assets')
ylabel('Lifetime Utility')
legend('High Income','Low Income') 


%% BASE CASE MODEL SIMULATION

% Initialize simulation
isim = 2*ones(nrep,1);                  % all agents employed,
ssim = zeros(nrep,1);                   % no assets
rng('default')

% Initialize history arrays
isimhist = zeros(nrep,nper);   % income state
ssimhist = zeros(nrep,nper);   % beginning assets
xsimhist = zeros(nrep,nper);   % change in assets

% Simulate model
for t=1:nper
  % store states in history arrays
  isimhist(:,t) = isim;
  ssimhist(:,t) = ssim;
  if t<nper
    % new employment state
    isimnew = discrand(nrep,w);
    % new ending assets
    ssimnew  = zeros(nrep,1);
    for i=1:ny
      ind = find(isim==i);
      if ~isempty(ind)
        ssimnew(ind) = interp1(s,send(:,i),ssim(ind));
      end
    end
    xsimhist(:,t) = ssimnew-ssim;
    isim = isimnew;
    ssim = ssimnew;
  end
end

% Plot simulated state path
figure
hold on
T = 30;
plot(0:T,ssimhist(1:3,1:T+1),0:T,mean(ssimhist(:,1:T+1)),'k')
ytickformat('%.1f')
title('Simulated and Expected Beginning Assets')
xlabel('Period')
ylabel('Assets')

% Retain post-warmup period observations 
isimhist = isimhist(:,nwrm:nper); isimhist = isimhist(:);
ssimhist = ssimhist(:,nwrm:nper); ssimhist = ssimhist(:);
xsimhist = xsimhist(:,nwrm:nper); xsimhist = xsimhist(:);

% Compute ergodic mean beginning assets, change in assets and visit probabilities, per income state
savg = zeros(ny,1);
xavg = zeros(ny,1);
prop = zeros(ny,1);
for i=1:ny
  savg(i) = mean(ssimhist(isimhist==i));
  xavg(i) = mean(xsimhist(isimhist==i));
  prop(i) = mean(isimhist==i);
end
xavgall = mean(xsimhist);
savgall = mean(ssimhist);

% Compute ergodic standard deviations of assets and asset changes, per income state
sstd = zeros(ny,1);
xstd = zeros(ny,1);
for i=1:ny
  sstd(i) = std(ssimhist(isimhist==i)+xsimhist(isimhist==i));
  xstd(i) = std(xsimhist(isimhist==i));
end
sstdall = std(ssimhist+xsimhist);
xstdall = std(xsimhist);

% Print ergodic mean assets and visit probabilities, per income state
fprintf('\nErgodic Mean Assets and Visit Probabilities, per Income State\n')
fprintf('                     High       Low\n')
fprintf('                   Income     Income       All\n')
fprintf('Ending Assets     %7.3f  %9.3f %9.3f\n',savg+xavg,savgall+xavgall)
fprintf('Beginning Assets  %7.3f  %9.3f %9.3f\n',savg,savgall)
fprintf('Change in Assets  %7.3f  %9.3f %9.3f\n',xavg,xavgall)
fprintf('Probability       %7.3f  %9.3f %9.3f\n',prop,sum(prop))

% Print ergodic standard deviation of assets, per income state
fprintf('\nErgodic Standard Devitations of Assets, per Income State\n')
fprintf('                     High       Low\n')
fprintf('                   Income     Income       All\n')
fprintf('Beginning Assets  %7.3f  %9.3f %9.3f\n',sstd,sstdall)
fprintf('Change in Assets  %7.3f  %9.3f %9.3f\n',xstd,xstdall)

end

%% SAVE FIGURES
printfigures(mfilename)


%% RESIDUAL FUNCTION
function [r,send] = resid(v,n,ny,util,s,S,Is,IS,y,w,tau,delta)
% Reshape and intialize
v  = reshape(v,n,ny);
v1 = zeros(n,ny);
v2 = zeros(n,ny);
send = zeros(n,ny);
% Evaluate RHS of Bellman equation
Evnext1 = v*w;
Evnext2 = Evnext1(:,ones(1,n))';
for i=1:ny
  % do not transact
  v1(:,i) = util(y(i)+Is) + delta*Evnext1;
  % transact
  C = y(i) + IS + S - S' - tau;
  U = util(C);
  U(C<=0) = -inf;
  [v2(:,i),isn] = max(U+delta*Evnext2,[],2);
  send(:,i) = s(isn);
end
r = v(:) - max(v1(:),v2(:));
for i=1:ny
  ind = find(v1(:,i)>v2(:,i));
  send(ind,i) = s(ind);
end
end