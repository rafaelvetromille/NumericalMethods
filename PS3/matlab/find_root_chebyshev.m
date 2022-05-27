% Criamos uma função que encontra as raízes do Polinômio de Chebychev de d+1 graus

function [roots, kgrid_roots] = find_root_chebyshev(d, parameters)
    
    kmin        = parameters.kgrid(1); 
    kmax        = parameters.kgrid(end);

    params_size = d + 1;
    roots       = zeros(params_size,1);
    kgrid_roots = zeros(params_size,1);
    
for i = 1:params_size
    
    roots(i)       = -cos((2*i-1)/(2*params_size) * pi);      % vide slide pag 19
    kgrid_roots(i) = ((1 + roots(i))/2)*(kmax-kmin) + kmin;   % retornando ao valor normal!

end

end