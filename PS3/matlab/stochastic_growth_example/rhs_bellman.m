function y = rhs_bellman(c,s,parameters,a,fspace)

beta  = parameters.beta;
sigma = parameters.sigma;

u  = utility(c,sigma);

Ev = expected_value(c,s,parameters,a,fspace);

y  = u+beta*Ev;

