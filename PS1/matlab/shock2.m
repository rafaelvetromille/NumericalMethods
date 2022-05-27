%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: th = shock(th0,Pi,T)
%
% This function simulates a time series of length T for a discrete markov
% process with transition matrix Pi and initial state th0.
%
% INPUTS
%   th0:    Initial state of the process at date 1
%   Pi:     Transition matrix for the discrete markov process. Rows are
%           current state and columns next period state
%   T:      Number of time periods to be simulated
%
% OUTPUTS
%   th:     Column vector with simulated time series. Contains the index of
%           the current state of the process
%   
% Author: Jan Hannes Lang
% Date: 28.04.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = shock2(th0, eps, Pi,T)

idx    = ones(T,1);
idx(1) = th0;
cum   = cumsum(Pi')';

for i = 2:T
    x     = find(normcdf(eps(i), 0, 0.007) <= cum(idx(i-1),:));
    idx(i) = x(1);
end
end
