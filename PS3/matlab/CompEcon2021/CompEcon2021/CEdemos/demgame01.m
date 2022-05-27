%% DEMGAME01 Production Capacity Game Model
%
% Two firms competing in a duopolistic product market must decide how much
% to invest in productive capacity.

%    States
%        s       capital stock at beginning of period
%    Actions
%        x       investent this period
%    Parameters
%        alpha   parameters
%        beta    parameters
%        gamma   parameters
%        psi     depreciation rate
%        delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';                     % basis function type
n     = [8 8];                          % degree of approximation
smin  = [0.5 0.5];                      % minimum state
smax  = [1.5 1.5];                      % maximum state
nper  = 30;                             % number of periods simulated
nrep  = 2;                              % number of replications

% Base case model parameters
alpha = [8.0 4.0];                      % parameters
beta  = [1.8 0.2];                      % parameters
gamma = [0.4 3.0];                      % parameters
psi   = 0.1;                            % depreciation rate
delta = 0.9;                            % discount factor

% Model structure
clear model
model.func = @func;                     % model functions
model.discount = delta;                 % discount factor
model.params = {alpha beta gamma psi};  % other parameters
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i

% Approximation structure
basis = fundefn(basistype,n,smin,smax);% basis functions

% Check model derivatives
gamecheck(model,(smax+smin)/2,zeros(1,2));


%% BASE CASE MODEL SOLUTION

% Solve collocation equations
optset('gamesolve','nres',4)
[s,v,c,sr,vr,xr,resid] = gamesolve(model,basis);

% Reshape output for plotting
nr = n*4;
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

% Plot optional investment policies
for ip=1:2
  figure
  mesh(sr1,sr2,xr(:,:,ip));
  xtickformat('%.1f')
  ytickformat('%.1f')
  ztickformat('%.2f')
  title(['Optimal Investment Policy: Player ' num2str(ip)])
  xlabel('$S_1$')
  ylabel('$S_2$')
  zlabel('Investment')
  zlim([0 .25])
end

% Plot value functions
for ip=1:2
  figure
  mesh(sr1,sr2,vr(:,:,ip));
  xtickformat('%.1f')
  ytickformat('%.1f')
  ztickformat('%.1f')
  title(['Value Function: Player ' num2str(ip)])
  xlabel('$S_1$')
  ylabel('$S_2$')
  zlabel('Value')
end

% Plot own shadow price functions
for ip=1:2
  figure
  mesh(sr1,sr2,pr(:,:,ip));
  xtickformat('%.1f')
  ytickformat('%.1f')
  ztickformat('%.1f')
  title(['Own Shadow Price: Player ' num2str(ip)])
  xlabel('$S_1$')
  ylabel('$S_2$')
  zlabel('Price')
end

% Plot residuals
for ip=1:2
  figure
  mesh(sr1,sr2,resid(:,:,ip));
  xtickformat('%.1f')
  ytickformat('%.1f')
  title(['Approximation Residual: Player ' num2str(ip)])
  xlabel('$S_1$')
  ylabel('$S_2$')
  zlabel('Residual')
end


%% BASE CASE MODEL SIMULATION

% Simulate model
ss    = cell(2,1);
ss{1} = sr1(:,1);
ss{2} = sr2(1,:)';
sinit = [smin(1)*ones(nrep,1) smax(2)*ones(nrep,1)];
esim = zeros(nrep,nper+1,2);
[spath,xpath] = gamesimul(model,sinit,nper,ss,xr,esim);
smean = squeeze(spath(1,:,:))';
xmean = squeeze(xpath(1,:,:))';

% Plot capital stock paths
figure
plot((0:nper)',smean)
ytickformat('%.1f')
legend('Player 1','Player 2')
title('Simulated Capital Stock')
xlabel('Period')
ylabel('Capital Stock')

% Plot investment paths
figure
plot((0:nper)',xmean)
ytickformat('%.2f')
legend('Player 1','Player 2')
title('Simulated Investment')
xlabel('Period')
ylabel('Investment')


%% SAVE FIGURES
printfigures(mfilename)


%% GAMESOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,i,s,x,~,alpha,beta,gamma,psi)
ns = size(s,1);
if ns>1
  x = squeeze(x);
end
switch flag
  case 'b'
    xl = zeros(ns,1);
    xu = inf*ones(ns,1);
    out1=xl; out2=xu; out3=[];
  case 'f'
    c    = beta(1) + beta(2)./s;
    prof = ((alpha(1)-2*c(:,i)+c(:,3-i)).^2)/(9*alpha(2));
    f    = prof-(gamma(1)*x(:,i)+0.5*gamma(2)*x(:,i).^2);
    fx   = -(gamma(1)+gamma(2)*x(:,i));
    fxx  = zeros(ns,1)-gamma(2);
    out1=f; out2=fx; out3=fxx;
  case 'g'
    gx  = zeros(ns,2,1);
    gxx = zeros(ns,2,1,1);
    g   = (1-psi)*s + x;
    gx(:,i) = ones(ns,1);
    out1=g; out2=gx; out3=gxx;
end
end