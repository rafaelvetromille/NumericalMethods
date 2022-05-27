%% QNWNORMUNIF
%
%  Generates a discrete n-valued uniformly-distributed approximation to a
%  univatiate normal random variable with mean mu and variance xvar.
%
%  Usage
%    [x,w] = qnwnormunif(n,mu,xvar)
%  Input
%    n    : degree of discretization
%    mu   : mean (default: 0)
%    xvar : variance (default: 1)
%  Output
%    x    : n.1 discrete mass points
%    w    : n.1 associated probabilities (each 1/n)
%  Note
%    This is not a Guassian discretization scheme.  Moments beyond the
%    second only approximately replicated by discrete approximant.
%  Note
%    To compute Ef(X) when f is real-valued on Re and X is univariate
%    Normal(mu,var), write a vectorized Matlab function f that returns an
%    m.1 vector when passed an m.1 vector, and execute
%    [x,w]=qnwnormunif(m,mu,var) and set Ef=w'*f(x) or Ef=mean(f(x)).

%  Copyright(c) 2010-2021
%   Mario J. Miranda - miranda.4@osu.edu

function [x,w] = qnwnormunif(n,mu,xvar)
if nargin<2, mu   = 1; end
if nargin<3, xvar = 1; end
if n==1||xvar==0
  x = mu*ones(n,1);
else
  c = nodeunif(n+1,0,1);
  x = icdf('Normal',(c(1:n)+c(2:n+1))/2,0,1);
  x = mu+sqrt((xvar/var(x,1)))*x;
end
w = ones(n,1)/n;