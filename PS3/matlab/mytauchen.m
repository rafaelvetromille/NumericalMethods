% FUNCTION: [s, Pi] = mytauchen(const,rho,sig,N)
%
% This function discretizes a continuous AR(1) process by using the method
% proposed by Tauchen (1986). The AR(1) process takes the following form:
% y(t) = const + rho*y(t-1) + eps(t), where eps ~ N(0,sig^2)
%
% INPUTS
%   const:  scalar, intercept term of the AR(1) process
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%   N:      scalar, number of grid points for the discretized process
%
% OUTPUTS
%   s:      column vector of size Nx1, contains all possible states in ascending order
%   Pi:     matrix of size NxN, contains the transition proabilities. Rows
%           are current state and columns future state

function [theta, Pi] = mytauchen(const,m,rho,sig,N)

s       = zeros(N,1);
Pi      = eye(N,N);

s(N)    = const/(1-rho) + m * sqrt(sig^2/(1-rho^2)); 
s(1)    = const/(1-rho) - m * sqrt(sig^2/(1-rho^2)); 

theta(N) = exp(s(N)); % theta_h
theta(1) = exp(s(1)); % theta_l

step    = (s(N)-s(1))/(N-1);

for i=2:(N-1)
   s(i) = s(i-1) + step; 
   theta(i) = exp(s(i));
end

for j = 1:N
    for k = 1:N
        if k == 1
            Pi(j,k) = cdf_normal((s(1) - const - rho*s(j) + step/2) / sig);
        elseif k == N
            Pi(j,k) = 1 - cdf_normal((s(N) - const - rho*s(j) - step/2) / sig);
        else
            Pi(j,k) = cdf_normal((s(k) - const - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - const - rho*s(j) - step/2) / sig);
        end
    end
end

function c = cdf_normal(x)
    c = (1/2) * (1 + erf(x/sqrt(2)));
end

end