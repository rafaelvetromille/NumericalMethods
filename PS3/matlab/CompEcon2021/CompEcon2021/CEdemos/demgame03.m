%% DEMGAME03 Marketing Board Game Model
%
% The marketing boards of two countries that compete in a duopolistic 
% international agricultural commodity market must decide how much to 
% export and how much to store.

%    States
%        s       quantity available per country at beginning of period
%    Actions
%        x       stocks hend in inventory per country at end of period
%    Parameters
%        kappa   unit storage cost
%        gamma   inverse demand elasticity
%        xmax    storage capacity
%        mu      relative size of agents
%        yvol    yield volatility
%        delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb'; % basis function type
n     = [20 20];                        % degree of approximation
m     = [5 5];                          % number of shocks
nper  = 20;                             % number of periods simulated
nrep  = 1000;                           % number of replications

% Base case model parameters
kappa = 0.05;                           % unit storage cost
gamma = -0.5;                           % inverse demand elasticity
xmax  = [0.2 0.2];                      % storage capacity
mu    = [0.5 0.5];                      % relative size of agents
yvol  = 0.2;                            % yield volatility
delta = 0.95;                           % discount factor

% Shock distribution
yvar  = (yvol^2)*eye(2);                % covariance matrix of log shocks
[e,w] = qnwlogn(m,mu,yvar);             % log shocks and probabilities

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {kappa gamma xmax};      % other parameters
model.discount = delta;                 % discount factor
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i
model.e = e;                            % shocks
model.w = w;                            % probabilities

% Approximation structure
smin   = min(e);                        % minimum supply
smax   = max(e)+xmax;                   % maximum supply
basis  = fundefn(basistype,n,smin,smax);   % basis functions

% Check model derivatives
gamecheck(model,(smax+smin)/2,zeros(1,2));


%% BASE CASE MODEL SOLUTION

% Solve Bellman Equation
optset('gamesolve','nres',3)
[v,x,c,sr,vr,xr,resid] = gamesolve(model,basis);

% Reshape output for plotting
nr = n*3;
sr1 = reshape(sr(:,1),nr);
sr2 = reshape(sr(:,2),nr);
vr = reshape(vr,[nr 2]);
xr = reshape(xr,[nr 2]);
resid = reshape(resid,[nr 2]);

% Compute shadow prices
pr = zeros(prod(nr),2);
pr(:,1) = funeval(c(:,1),basis,sr,[1 0]);
pr(:,2) = funeval(c(:,2),basis,sr,[0 1]);
pr = reshape(pr,[nr 2]);

% Plot optimal policies
for ip=1:2
  figure
  mesh(sr1,sr2,xr(:,:,ip))
  xlim([sr1(1),sr1(end)])
  ylim([sr2(1),sr2(end)])
  zlim([0 0.12])
  ztickformat('%.2f')
  title(['Optimal Ending Stocks: Country ' num2str(ip)])
  xlabel('$S_1$');
  ylabel('$S_2$')
  zlabel('Ending Stocks')
end

% Plot value functions
for ip=1:2
  figure
  mesh(sr1,sr2,vr(:,:,ip))
  xlim([sr1(1),sr1(end)])
  ylim([sr2(1),sr2(end)])
  ztickformat('%.1f')
  title(['Value Function: Country ' num2str(ip)])
  xlabel('$S_1$');
  ylabel('$S_2$')
  zlabel('Value')
end

% Plot own shadow price functions
for ip=1:2
  figure
  mesh(sr1,sr2,pr(:,:,ip))
  xlim([sr1(1),sr1(end)])
  ylim([sr2(1),sr2(end)])
  title(['Own Shadow Price: Country ' num2str(ip)])
  xlabel('$S_1$');
  ylabel('$S_2$')
  zlabel('Price')
end

% Plot market price function
figure
sg = gridmake(sr);
xx = reshape(xr,size(sg));
pp = sum(sg-xx,2).^gamma;
pp = reshape(pp,nr);
mesh(sr1,sr2,pp)
xlim([sr1(1),sr1(end)])
ylim([sr2(1),sr2(end)])
zlim([0.6 1.4])
ztickformat('%.1f')
title('Market Price')
xlabel('$S_1$');
ylabel('$S_2$')
zlabel('Market Price')

% Plot residuals
for ip=1:2
  figure
  mesh(sr1,sr2,resid(:,:,ip))
  xlim([sr1(1),sr1(end)])
  ylim([sr2(1),sr2(end)])
  ztickformat('%.1f') 
  title(['Residual: Country ' num2str(ip)])
  xlabel('$S_1$');
  ylabel('$S_2$')
  zlabel('Residual')
end


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(mu,yvar,nrep,nper+1);

% Initialize simulation
ss{1} = sr1(:,1);
ss{2} = sr2(1,:)';
sinit = [smin(1)*ones(nrep,1) smax(2)*ones(nrep,1)];
xr = reshape(xr,[nr 1 2]);

% Simulate model
[spath,xpath] = gamesimul(model,sinit,nper,ss,xr,esim);
smean = squeeze(mean(spath,1))';
xmean = squeeze(mean(xpath,1))';

% Plot expected supply paths
for ip=1:2
  figure
  hold on
  plot(0:nper,spath(1:3,:,ip))
  plot(0:nper,mean(spath(:,:,ip)),'k')
  ytickformat('%.1f')
  title(['Simulated and Expected Supply: Country ' num2str(ip)])
  xlabel('Period')
  ylabel('Supply')
end

% Plot expected ending stocks paths
for ip=1:2
  figure
  hold on
  plot(0:nper,xpath(1:3,:,ip))
  plot(0:nper,mean(xpath(:,:,ip)),'k')
  ytickformat('%.2f')
  title(['Simulated and Expected Stocks: Country ' num2str(ip)])
  xlabel('Period')
  ylabel('Ending Stocks')
end


%% SAVE FIGURES
printfigures(mfilename)


%% GAMESOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,i,s,x,e,kappa,gamma,xmax)
ns = size(s,1);
if ns>1
  x = squeeze(x);
end
switch flag
  case 'b'
    xl = zeros(ns,1);
    xu = xmax(i)*ones(ns,1);
    out1=xl; out2=xu; out3=[];
  case 'f'
    q   = s-x;
    qtot = sum(q,2);
    p   = qtot.^gamma;
    px  = -gamma*qtot.^(gamma-1);
    pxx = (gamma-1)*gamma*qtot.^(gamma-2);
    f   = p.*q(:,i) - kappa*x(:,i);
    fx  = px.*q(:,i) - p - kappa;
    fxx = pxx.*q(:,i) - 2*px;
    out1=f; out2=fx; out3=fxx;
  case 'g'
    gx  = zeros(ns,2,1);
    gxx = zeros(ns,2,1,1);
    g   = x + e;
    gx(:,i) = ones(ns,1);
    out1=g; out2=gx; out3=gxx;
end
end