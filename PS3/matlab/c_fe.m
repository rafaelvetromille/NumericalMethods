function cons = c_fe(a, k, n_elements, parameters)

sum   = 0;

for ik = 1:n_elements
    sum = sum + a(ik)*psi_function(ik, k, n_elements, parameters);
end

cons = sum;

end
