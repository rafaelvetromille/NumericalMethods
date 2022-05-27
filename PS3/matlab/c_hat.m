function cons = c_hat(gamma, k, d, parameters)

% kmax e kmin
kmin     = parameters.kgrid(1);
kmax     = parameters.kgrid(end);

% transladar o grid de k
k_trans = 2*(k - kmin)/(kmax - kmin) - 1;   
order   = d + 1;
sum     = 0;

for i = 1:order
    sum = sum + gamma(i)*chebyshev(i-1, k_trans);
end

cons = sum;

end
