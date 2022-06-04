function cons = c_fe(a, k, n_elements, parameters)

sum   = 0;

for in = 1:n_elements
    sum = sum + a(in)*psi_function(in, k, n_elements, parameters);
end

cons = sum;

end
