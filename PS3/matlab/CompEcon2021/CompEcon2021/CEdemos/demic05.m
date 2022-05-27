%% DEMIC05 Cash Management Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION
  
% Base case model parameters
mu    = 0;
sigma = 0.2;
rho   = 0.04;
r     = 0.05;
c     = 0.75;
C     = 0.25;
f     = 1;
F     = 0.5;

% Model structure
clear model
model.func   = @func;
model.params = {mu sigma rho r c C};
model.xindex = [1 1;2 1];
model.F      = [f;F];


%% BASE CASE MODEL SOLUTION

% Initial values
x = [0 1;2 1];

% Solve collocation equation
n = 10;
[cv,basis,x] = icsolve(model,x,n);

sstar  = x(:);
vstar  = funeval(cv,basis,sstar);
vstar1 = funeval(cv,basis,sstar,1);

% Value function and derivatives
S  = nodeunif(101,0,x(2,1));
V  = funeval(cv,basis,S);
Vs = funeval(cv,basis,S,1);

% Plot value function
figure
hold on
plot(S,V)
plotvdash(sstar(1),vstar(1))
plotbullet(sstar(1),vstar(1))
plottext(sstar(1)+0.01,[],'$s_1^*$')
plotvdash(sstar(2),vstar(2))
plotbullet(sstar(2),vstar(2))
plottext(sstar(2)+0.01,[],'$s_2^*$')
xtickformat('%.1f')
title('Value Function')
xlabel('$S$')
ylabel('$V$')

% Plot first derivative of value function
figure 
hold on
plot(S,Vs)
plothdash([],ones(2,1)*(r/rho+[c -C]))
plotbullet(sstar,vstar1)
xtickformat('%.1f')
title('First Derivative of Value Function')
xlabel('$S$')
ylabel('$V''$')

% Display selected values
fprintf('\n') 
fprintf('Switching Points\n') 
fprintf('           x     V(x)    V''(x)\n')
fprintf('      %6.3f   %6.3f    %5.3f\n',x(1),vstar(1),vstar1(1))
fprintf('      %6.3f   %6.3f    %5.3f\n',x(3),vstar(3),vstar1(3))
fprintf('      %6.3f   %6.3f    %5.3f\n',x(4),vstar(4),vstar1(4))
fprintf('      %6.3f   %6.3f    %5.3f\n',x(2),vstar(2),vstar1(2))


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
function out = func(flag,s,mu,sigma,rho,r,c,C)
switch flag
  case 'f'        % reward
    out = zeros(size(s,1),1);
  case 'mu'       % state drift
    out = mu+zeros(size(s,1),1);
  case 'sigma'    % state diffusion
    out = sigma+zeros(size(s,1),1);
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = [(-c-r/rho)*(s(2)-s(1)) (c+r/rho) (-c-r/rho) 0];
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    out = [(r/rho-C)*(s(1)-s(2)) (r/rho-C) (C-r/rho) 0];
end
end