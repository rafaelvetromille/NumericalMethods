%% DISCMOMENTS
%
%  Generates mean and covariance of discrete multivariate distribution.
%
%  Usage
%    [xavg,xvar] = discmoments(x,w)
%  Let
%    d = number of random variates
%    n = number of probability mass points
%  Input
%    x    : n.d values for d distinct variates
%    w    : n.1 discrete probabilities
%  Output
%    xavg : 1.d mean vector
%    xvar : d.d variance-covariance matrix

%  Copyright(c) 2014-2021
%   Mario J. Miranda - miranda.4@osu.edu

function [xavg,xvar] = discmoments(x,w)
xavg = w'*x;
d = size(x,2);
if d==1
  xvar = w'*((x-xavg).^2);
else
  xvar = zeros(d,d);
  for i=1:d
    for j=1:d
      xvar(i,j) = w'*((x(:,i)-xavg(i)).*(x(:,j)-xavg(j)));
    end
  end
end