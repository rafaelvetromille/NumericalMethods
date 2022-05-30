function residual_vector = build_system_fe(x, nl, parameters)

    nz           = length(parameters.zgrid);    
    matr_aux     = zeros(nz,nl);
    
    kmax         = parameters.kgrid(end);
    kmin         = parameters.kgrid(1);
    kgrid_aux    = linspace(kmin, kmax, nl);
    
for iz = 1:nz
    for il = 1:nl
        matr_aux(iz,il) = residual_function_fe(x, kgrid_aux(il), nl, iz, parameters);
    end
end

    residual_vector = reshape(matr_aux, [], 1);

end