%% DEMDDP08 Deterministic Cow Replacement Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Base case model parameters
delta = 0.9;                            % discount factor
cost  = 500;                            % replacement cost
price = 150;                            % milk price

% State space
S = (1:10)';                            % lactation states
n = length(S);                          % number of states

% Action space (keep='K', replace='R')
X = ['K';'R'];                          % keep or replace
m = length(X);                          % number of actions

% Reward function
y = (-0.2*S.^2+2*S+8);                  % yield per lactation
f = [price*y price*y-cost];             % net revenue by action
f(10,1) = -inf;                         % force replace at lactation 10

% State transition function
g = zeros(n,m);
for i=1:n
  g(i,1) = min(i+1,n);                % increase lactation number by 1, if keep
  g(i,2) = 1;                         % lactation number reverts to 1, if replace
end

% Model structure
clear model
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% BASE CASE MODEL SOLUTION

% Solve Bellman equation
[v,x] = ddpsolve(model);

% Plot optimal policy
figure
bar(S,x)
yticks([1 2])
set(gca,'YTickLabel',{'Keep','Replace'})

% Plot value function
figure
plot(S,v)
ytickformat('%.2f')
title('Value Function')
xlabel('Age')
ylabel('Value of Operation')


%% SAVE FIGURES
printfigures(mfilename)