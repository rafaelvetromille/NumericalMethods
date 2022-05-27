%% DEMDDP05 Bioeconomic Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Base case model parameters
T     = 10;                             % foraging periods
emax  =  8;                             % energy capacity
e     = [2  4  4];                      % energy from foraging
p     = [1.0 0.7 0.8];                  % predation survival probabilities
q     = [0.5 0.8 0.7];                  % foraging success probabilities

% State space
S = (0:emax)';                          % energy levels
n = length(S);                          % number of states

% Action space
X = (1:3)';                             % vector of actions
m = length(X);                          % number of actions

% Reward function
f = zeros(n,m);

% State transition probability matrix
P = zeros(m,n,n);
for k=1:m
  P(k,1,1) = 1;
  for i=2:n
    % does not survive predation
    snext = 0;           j=getindex(snext,S); P(k,i,j) = P(k,i,j) + (1-p(k));
    % survives predation, finds food
    snext = S(i)-1+e(k); j=getindex(snext,S); P(k,i,j) = P(k,i,j) + p(k)*q(k);
    % survives predation, finds no food
    snext = S(i)-1;      j=getindex(snext,S); P(k,i,j) = P(k,i,j) + p(k)*(1-q(k));
  end
end

% Terminal value function
vterm = ones(n,1);                      % terminal value: survive
vterm(1) = 0;                           % terminal value: death

% Model structure
clear model
model.reward     = f;
model.transprob  = P;
model.horizon    = T;
model.discount   = 1;
model.vterm      = vterm;


%% BASE CASE MODEL SOLUTION

% Solve Bellman equation
[v,x] = ddpsolve(model);

% Plot survival probabilities, Period 0
figure
bar(S,v(:,1),1)
axis([-.5 emax+0.5 0 1]);
ytickformat('%.1f')
title('Survival Probability (Period 0)')
xlabel('Stock of Energy')
ylabel('Probability')

% Plot survival probabilities, Period 5
figure
bar(S,v(:,6),1)
axis([-.5 emax+0.5 0 1])
ytickformat('%.1f')
title('Survival Probability (Period 5)')
xlabel('Stock of Energy')
ylabel('Probability')


%% SAVE FIGURES
printfigures(mfilename)