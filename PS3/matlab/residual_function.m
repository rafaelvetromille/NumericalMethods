% Resiual function 

function r_calc = residual_function(gamma, k, d, j, parameters)

nz           = length(parameters.zgrid);

co           = c_hat(gamma(j,:), k, d, parameters);
g            = find_g(k, parameters.zgrid(j), co, parameters);

aux1         = zeros(nz,1);
aux2         = zeros(nz,1);

for s = 1:nz

    c_prime = c_hat(gamma(s,:), g, d, parameters);
    
    aux1(s) = (1 - parameters.delta + parameters.alpha*parameters.zgrid(s)*g^(parameters.alpha-1));
    aux2(s) = (c_prime/co)^(-parameters.mu);

end

aux       = aux1 .* aux2;
r_calc    = parameters.beta * parameters.P(j,:)*aux - 1;

% r_calc = beta * (c_prime/co)^(-mu) * (1 - delta + alpha*z*k_prime^(alpha-1)) - 1  % caso discreto
end