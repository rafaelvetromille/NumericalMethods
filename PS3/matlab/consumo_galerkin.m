function cons = consumo_galerkin(x,k, n_k)



soma = 0;

for j = 1:n_k

    soma = soma + x(j)*psi_i(j,k, n_k);

end

cons = soma;

end

