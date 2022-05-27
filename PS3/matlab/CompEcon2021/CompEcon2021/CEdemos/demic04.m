%% DEMIC04 Capacity Choice Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION
  
% Base case model parameters
P     = 2;
C     = 1;
delta = 0.5;
rho   = 0.1;

% Model structure
clear model
model.func   = @func;
model.params = {P C delta rho};
model.xindex = [2 0;0 0];
model.F      = [0;0];


%% BASE CASE MODEL SOLUTION

% Initial values
xmax = 20;
x = [0 0;xmax 0];

% Solve collocation equation
n = 25;
[cv,basis,x] = icsolve(model,x,n);
Kstar = x(1,1);

% Value function and derivatives
K = nodeunif(101,0,xmax);
V = funeval(cv,basis,K);
Vstar = funeval(cv,basis,Kstar);
V(K<Kstar) = Vstar-(Kstar-K(K<Kstar))*C;

% Plot value function
figure
hold on
plot(K,V)
plotvdash(Kstar,Vstar)
plotbullet(Kstar,Vstar)
plottext(Kstar+0.01,[],'$K^*$')
title('Value Function')
xlabel('K')
ylabel('V')

% Display selected values
fprintf('\n') 
fprintf('Switching Points\n') 
fprintf('         K*    V(K*)\n')
fprintf('      %5.3f   %5.3f\n',Kstar,Vstar)


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
function out = func(flag,s,P,C,delta,rho)
switch flag
  case 'f'        % reward
    out = P*log(s+1);
  case 'mu'       % state drift
    out = -delta*s;
  case 'sigma'    % state diffusion
    out = [];
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = [C*(s(1)-s(2)) C -C 0];
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    out = zeros(1,4);
end
end