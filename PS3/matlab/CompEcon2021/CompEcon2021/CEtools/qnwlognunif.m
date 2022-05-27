%% QNWLOGNUNIF
%
%  Generates a discrete n-valued uniformly-distributed approximation to a
%  univatiate lognormal random variable with mean mu and log variance
%  logvar.
%
%  Usage
%    [x,w] = qnwlognunif(n,mu,logvar)
%  Input
%    n      : degree of discretization
%    mu     : mean of x (default: 1)
%    logvar : log variance of x (default: 1)
%  Output
%    x   : n discrete mass points
%    w   : n associated probabilities (each 1/n)
%  Note
%    This is not a Guassian discretization scheme.  Moments beyond the
%    second only approximately replicated by discrete approximant.
%  Note
%    To compute Ef(X) when f is real-valued on Re and X is univariate
%    lognormal random variable with mean mu and log variance logvar, write
%    a vectorized Matlab function f that returns an m.1 vector when passed
%    an m.1 vector of values, execute [x,w]=qnwlognunif(n,mu,logvar), and
%    set Ef=w'*f(x) or Ef=mean(f(x)).

%  Copyright(c) 2021
%   Mario J. Miranda - miranda.4@osu.edu

function [x,w] = qnwlognunif(n,mu,logvar)
if nargin<2, mu = 1; end
if nargin<3, logvar = 1; end
if n==1||nargin<3||logvar==0
  x = mu*ones(n,1);
else
  x = qnwnormunif(n,0,logvar);
  x = exp(x);
  x = x*mu/mean(x);
end
w = ones(n,1)/n;