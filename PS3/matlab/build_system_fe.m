function residual_vector = build_system_fe(x, n_elements, parameters)

    nz           = length(parameters.zgrid);    
    matr_aux     = zeros(nz,n_elements);
    
    kmax         = parameters.kgrid(end);
    kmin         = parameters.kgrid(1);
     
%     index        = linspace(0, 0, n_elements); 
%     
%     for i = 2:(n_elements-1)
%        d            = floor(length(parameters.kgrid)/n_elements);
%        index(i)     = parameters.kgrid(1 + d*i);  
%     end
%     
%     index(1)     = kmin; 
%     index(end)   = kmax; 

    %kgrid_aux   = parameters.kgrid(index);
    kgrid_aux   = linspace(kmin, kmax, n_elements);
    
for iz = 1:nz
    for il = 1:n_elements
        matr_aux(iz,il) = residual_function_fe(x, kgrid_aux(il), n_elements, iz, parameters);
    end
end

    residual_vector = reshape(matr_aux, [], 1);

end