%% DEMODE05 Commodity Storage Model

%  Solve
%    s' = -p^(-eta)
%    p' = r*p+k
%  where
%    s: stocks
%    p: price
%    x = [s;p]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Velocity Function
r   = 0.1;     % interest rate
k   = 0.5;     % unit cost of storage
eta = 5;       % demand elasticity
f = @(x) [-x(2,:).^(-eta);r*x(2,:)+k];

% Time Horizon
T = 1;

% Boundary Conditions
s0  = 1;        % initial stocks
sT  = 0;        % terminal stocks
bx = [1;1];     % boundary variables
bt = [0;T];     % boundary times
bv = [s0;sT];   % boundary values


%% SOLVE ODE USING COLLOCATION METHOD (ODECOL)

% Solve ODE
n   = 15;        % number of basis functions
c = zeros(n,2); c(1,:) = 1;
[x,t,res] = odecol(f,bv,T,n,bt,bx,[],c);


%% PLOT SOLUTION

% Plot Solution
figure
plot(t,x)
ylim([0 1.5])
xtickformat('%.1f')
ytickformat('%.1f')
title('Commodity Storage Model')
xlabel('Time')
legend('Stocks','Price')

% Plot Residuals
figure
hold on
plot(t,res)
plothdash([],0)
xtickformat('%.1f')
title('Commodity Storage Model')
xlabel('Time')
ylabel('Residual')
legend('Stocks','Price')


%% SAVE FIGURES
printfigures(mfilename)