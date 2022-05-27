%% DEMOPT07 Bound-Constrained Optimization via Sequential LCP

% Preliminary tasks
deminit(mfilename)

% Generate problem test data
z = randn(2,2)-1;
a = min(z,[],2);
b = max(z,[],2);
x = rand(2,1);

% Set convergence parameters
maxit = 200;
eps = 1e-10;

% Perform sequential LCP
fprintf('\n')
for it=1:maxit
  xold = x;
  [f,d,s] = func(x);
  x = lcpsolve(s,d-s*xold,a,b,xold);
  change = norm(x-xold,inf);
  fprintf('%3i %10.2e\n',it,change)
  if change<eps, break, end
end

% Print results
fprintf('\n')
if it>=maxit
  warning('Sequential lcp failed in demomax');
else
  fprintf('Perform Bounded Maximization with Random Data\n');
  fprintf('       a         x         b        f''\n');
  fprintf('%8.2f  %8.2f  %8.2f  %8.2f\n',[a x b d]')
end


function [y,d,s] = func(x)
d = zeros(2,1);
s = zeros(2,2);
y = exp(-x(1)) + x(1)*x(2)^2;
d(1) = -exp(-x(1))+x(2)^2;
d(2) = 2*x(1)*x(2);
s(1,1) = exp(x(1));
s(1,2) = 2*x(2);
s(2,1) = 2*x(2);
s(2,2) = 2*x(1);
end