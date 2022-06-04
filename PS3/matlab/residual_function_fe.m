% Residual function 

function r_calc = residual_function_fe(a, k, n_elements, iz, parameters)

nz        = length(parameters.zgrid);

C0        = c_fe(a(iz,:), k, n_elements, parameters);
K1        = parameters.zgrid(iz)*(k^(parameters.alpha)) + (1-parameters.delta)*k - C0;

aux1      = zeros(nz,1);
aux2      = zeros(nz,1);

for s = 1:nz
     
    C1      = c_fe(a(s,:), K1, n_elements, parameters);
    
    aux1(s) = (1 - parameters.delta + parameters.alpha*parameters.zgrid(s)*K1^(parameters.alpha-1));
    aux2(s) = (C1/C0)^(-parameters.mu);
      
end

aux       = aux1 .* aux2;
r_calc    = parameters.beta * parameters.P(iz,:)*aux - 1;

end