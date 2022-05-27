% Criamos uma função que, dado K, nos entrega K'

function k_prime = find_g(k, z, c_func, parameters)

    k_prime = z*(k^parameters.alpha) + (1-parameters.delta)*k - c_func;
   
end