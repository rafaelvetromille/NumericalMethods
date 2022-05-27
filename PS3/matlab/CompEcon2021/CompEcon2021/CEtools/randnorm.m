%% RANDNORM
%
%  Uses Matlab routine randn to generate n.m instances of pseudo-random
%  multivariate (d-dimensional) normal random vector with 1.d mean vector
%  mu and d.d positive definite variance-covariance matrix var.
%
%  Usage
%    x = randnorm(mu,var,n,m)
%  Let
%    d = dimension of random vector
%  Input
%    mu  : 1.d mean vector
%    var : d.d positive definite variance matrix
%    n   : "rows" of pseudo-random vector generated
%    m   : "columns" of pseudo-random vector generated (default: 1)
%  Output
%    x   : n.m.d pseudo-random variates
%  Note
%    To compute Monte Carlo approximation of Ef(X) when f is real-valued on
%    Re^d and X is Normal(mu,var) on R^d, write a vectorized Matlab
%    function f that returns an m.1 vector when passed an m.d matrix, and
%    execute x=randnorm(mu,var,N); Ef=mean(f(x)) where N is a large number.

%  Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu

function x = randnorm(mu,var,n,m)
d = length(mu);
if nargin<4, m=1; end
x = mu(ones(n*m,1),:);
if max(var(:))>0
  z = randn(n*m,d);
  x = x + z*chol(var);
end
x = squeeze(reshape(x,n,m,d));
end