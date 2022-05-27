function residual_vector  = build_system(gamma, d, parameters)

    nz           = length(parameters.zgrid);

    params_size  = d + 1;
    matr_aux     = zeros(nz, params_size);
    [~, k_roots] = find_root_chebyshev(d, parameters);
    
for iz = 1:nz
    for ip = 1:params_size
        matr_aux(iz,ip) = residual_function(gamma, k_roots(ip), d, iz, parameters);
    end
end
    
    % Necessário empilhar para funcionar:
    residual_vector = reshape(matr_aux, [], 1);

end