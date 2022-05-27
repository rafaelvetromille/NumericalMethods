%% DEMQUA05 Expected Utility Insurance Model
%
% Willingness to pay for complete insurance against income shocks.

% Preliminary tasks
deminit(mfilename)

% Parameters
ybar  = 1.0;    % Expected income
var   = 0.1;    % Variance of log income
alpha = 2;      % Constant coefficient or relative risk aversion

n = 20;
[y,w] = qnwlogn(n,ybar,var);
expectedutility = -w'*exp(-alpha*y);
certainutility  = -exp(-alpha*ybar);
ystar = -log(-expectedutility)/alpha;
wtp = ybar-ystar;

fprintf('\n')
fprintf('Willingness to Pay for Complete Insurance\n')
fprintf('Expected Income             %6.3f\n',ybar)
fprintf('Certainty Equivalent Income %6.3f\n',ystar)
fprintf('Willingness to Pay          %6.3f\n',wtp)