function valor = psi_i(i,k, n_k)

%%%%%%%%%% DEFININDO PARÂMETROS %%%%%%%

alpha = 1/3;
beta = 0.987;
delta = 0.012;

%%%%%%% TAMANHO DO GRID DO CAPITAL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% CRIAR O GRID DO CAPITAL %%%%%%%%%%%%

k_ss = (1/beta + delta - 1)*(1/alpha);
k_ss = k_ss^(1/(alpha - 1));

grid_k = linspace(0.75*k_ss, 1.25*k_ss, n_k);

if i == 1
    
    if k >= grid_k(i) & k <= grid_k(i+1)
    
    valor = (grid_k(i+1) - k)/(grid_k(i+1) - grid_k(i));
    
    else
        
        valor = 0;
        
    end
    
elseif i == n_k
    
    if k <= grid_k(i) & k >= grid_k(i-1)
    
    valor = (k - grid_k(i-1))/(grid_k(i) - grid_k(i-1));
    
    else
        
        valor = 0;
        
    end
    
else
        
       if k >= grid_k(i-1) & k <= grid_k(i) 
           
            valor = (k - grid_k(i-1))/(grid_k(i) - grid_k(i-1));
            
       elseif k > grid_k(i) & k <= grid_k(i+1)
           valor = (grid_k(i+1) - k)/(grid_k(i+1) - grid_k(i));
           
       else 
           valor = 0;
           
       end
end



end

