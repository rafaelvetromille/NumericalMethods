function residual_vector = build_system_galerkin(a, nP, nI, P, zgrid, kgrid, alpha, beta, delta, mu, nZ)

%%% CONTADOR 
count           = 1;
residual_vector = zeros(nP*nI, 1)';

%%% 
kgrid_aux = linspace(kgrid(1), kgrid(end), nP);

for iz = 1:nZ
    
    for ip = 1:nP
        
        %%% DEFININDO GRID PARA AS INTEGRAIS
        if ip == 1
            
            grid_integral_lower = linspace(kgrid_aux(ip), kgrid_aux(ip+1), nI);
            
            aL                  = grid_integral_lower(1);
            bL                  = grid_integral_lower(nI);
            
        elseif ip == nP
            
            grid_integral_upper = linspace(kgrid_aux(ip-1), kgrid_aux(ip), nI);
            
            aU                  = grid_integral_upper(1);
            bU                  = grid_integral_upper(nI);
            
        else
            
            grid_integral_lower = linspace(kgrid_aux(ip), kgrid_aux(ip+1), nI);
            grid_integral_upper = linspace(kgrid_aux(ip-1), kgrid_aux(ip), nI);
            
            aL                  = grid_integral_lower(1);
            bL                  = grid_integral_lower(nI);
            
            aU                  = grid_integral_upper(1);
            bU                  = grid_integral_upper(nI);
            
        end
        
        %%%% CALCULAR AS INTEGRAIS:
        
        %%% PRIMEIRO CASO: EXTREMOS
        
        %%% LOWER BOUND
        if ip == 1
            
            aux0 = 0;
            
            for i = 1:nI
                
                %%% COMPUTAR O X, CONFORME JUDD (xi)
                xi    = cos((2*i-1)/(2*nI)*pi);
                
                %%% COMPUTAR O CAPITAL (K0)
                K0    = (1 + xi)*(bL - aL)/2 + aL;
                
                %%% CALCULAR O  CONSUMO (C0)
                C0    = consumo_galerkin(a(iz,:), K0, nP);
                
                %%% CALCULAR O KPRIME (K1)
                K1    = zgrid(iz)*(K0^(alpha)) + (1-delta)*K0 - C0;
                
                %%% CALCULAR A FUNÇÃO RESÍDUO
                aux11 = zeros(nZ,1); aux21 = zeros(nZ,1);
                
                for w = 1:nZ
                    C1       = consumo_galerkin(a(w,:), K1, nP);
                    aux11(w) = (1 - delta + alpha*zgrid(w)*K1^(alpha-1));
                    aux21(w) = (C1/C0)^(-mu);
                end
                
                aux    = aux11 .* aux21;
                r_calc = beta * P(iz,:)*aux - 1;
                
                %%% COMPUTAR INTEGRAL (SOMA DISCRETA)
                aux0   = aux0 + r_calc*psi_i(ip, K0, nP)*sqrt(1 - xi^2);
            end
            
            %%% ACUMULAR OS RESÍDUOS
            residual_vector(count) = pi*(bL - aL)*aux0/(2*nI);
            count = count + 1;
            
            %%% UPPER BOUND
        elseif ip == nP
            
            aux0 = 0;
            
            for i = 1:nI
                
                %%% COMPUTAR O X, CONFORME JUDD (xi)
                xi    = cos((2*i-1)/(2*nI)*pi);
                
                %%% COMPUTAR O CAPITAL (K0)
                K0    = (1 + xi)*(grid_integral_upper(nI) - grid_integral_upper(1))/2 + grid_integral_upper(1);
                
                %%% CALCULAR O  CONSUMO (C0)
                C0    = consumo_galerkin(a(iz,:), K0, nP);
                
                %%% CALCULAR O KPRIME (K1)
                K1    = zgrid(iz)*(K0^(alpha)) + (1-delta)*K0 - C0;
                
                %%% CALCULAR A FUNÇÃO RESÍDUO
                aux11 = zeros(nZ,1); aux21 = zeros(nZ,1);
                
                for w = 1:nZ
                    C1       = consumo_galerkin(a(w,:), K1, nP);
                    aux11(w) = (1 - delta + alpha*zgrid(w)*K1^(alpha-1));
                    aux21(w) = (C1/C0)^(-mu);
                end
                
                aux    = aux11 .* aux21;
                r_calc = beta * P(iz,:)*aux - 1;
                
                %%% COMPUTAR INTEGRAL (SOMA DISCRETA)
                aux0   = aux0 + r_calc*psi_i(ip, K0, nP)*sqrt(1 - xi^2);
                
            end
            
            %%% ACUMULAR OS RESÍDUOS
            residual_vector(count) = pi*(bU - aU)*aux0/(2*nI);
            count = count + 1;
            
            %%% SEGUNDO CASO: MEIO
            
        else
            
            %%% INTEGRAL 1
            aux1 = 0;
            
            for i = 1:nI
                
                %%% COMPUTAR O X, CONFORME JUDD (xi)
                xi    = cos((2*i - 1)/(2*nI)*pi);
                
                %%% CALCULAR O CAPITAL (K0)
                K0    = (1 + xi)*(bL - aL)/2 + aL;
                
                %%% CALCULAR O CONSUMO (C0)
                C0    = consumo_galerkin(a(iz,:), K0, nP);
                
                %%% CALCULAR O KPRIME (K1)
                K1    = zgrid(iz)*(K0^(alpha)) + (1-delta)*K0 - C0;
                
                %%% CALCULAR A FUNÇÃO RESÍDUO
                aux11 = zeros(nZ,1); aux21 = zeros(nZ,1);
                
                for w = 1:nZ
                    C1       = consumo_galerkin(a(w,:), K1, nP);
                    aux11(w) = (1 - delta + alpha*zgrid(w)*K1^(alpha-1));
                    aux21(w) = (C1/C0)^(-mu);
                end
                
                aux    = aux11 .* aux21;
                r_calc = beta * P(iz,:)*aux - 1;
                
                %%% COMPUTAR INTEGRAL (SOMA DISCRETA)
                aux1   = aux1 + r_calc*psi_i(ip, K0, nP)*sqrt(1 - xi^2);
                
            end
            
            %%% INTEGRAL 2
            aux2 = 0;
            
            for i = 1:nI
                
                %%% COMPUTAR O X, CONFORME JUDD (xi)
                xi = cos((2*i - 1)/(2*nI)*pi);
                
                %%% COMPUTAR O CAPITAL (K0)
                K0 = (1 + xi)*(bU - aU)/2 + aU;
                
                %%% CALCULAR O  CONSUMO (C0)
                C0 = consumo_galerkin(a(iz,:), K0, nP);
                
                %%% CALCULAR O KPRIME (K1)
                K1 = zgrid(iz)*(K0^(alpha)) + (1-delta)*K0 - C0;
                
                %%% CALCULAR A FUNÇÃO RESÍDUO
                aux11         = zeros(nZ,1); aux21         = zeros(nZ,1);
                
                for w = 1:nZ
                    C1       = consumo_galerkin(a(w,:), K1, nP);
                    aux11(w) = (1 - delta + alpha*zgrid(w)*K1^(alpha-1));
                    aux21(w) = (C1/C0)^(-mu);
                end
                
                aux    = aux11 .* aux21;
                r_calc = beta * P(iz,:)*aux - 1;
                
                %%% COMPUTAR INTEGRAL (SOMA DISCRETA)
                aux2   = aux2 + r_calc*psi_i(ip, K0, nP)*sqrt(1 - xi^2);
                
            end
            
            %%% ACUMULAR OS RESÍDUOS
            residual_vector(count) = pi*(bL - aL)*aux1/(2*nI) + pi*(bU - aU)*aux2/(2*nI);
            count    = count + 1;
            
        end
        
    end
    
end

end

