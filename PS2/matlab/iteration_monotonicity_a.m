function [Tv, idx] = iteration_monotonicity_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

iter  = 0;
error = 1;

nk    = length(kgrid);
nz    = length(zgrid);

Tv    = zeros(nk,nz);
idx   = zeros(nk,nz);

% Iteração da Função valor - Monotonicidade c/ Acelerador
% -------------------------------------------------------
tic;
while (error > tol && iter <= max_iter)
    
    if iter < 100 || rem(iter,10) == 0
        
        for iz = 1:nz
            
            z = zgrid(iz);
            c      = zeros(nk,1);
            H      = zeros(nk,1);
            
            % o piso (floor) vai sendo atualizado a cada iteração
            % começamos pelo caso igual ao anterior
            % ---------------------------------------------------
            floor = 1;
            
            for ik = 1:nk
                
                % Pegando o capital ik
                % -----------------------------------------------------------
                k      = kgrid(ik);
                
                % Esse loop define o vetor de consumo de forma a garantir que
                % nenhuma entrada será negativa: a cada elemento, tiramos o
                % máximo entre o consumo e zero
                
                % Note que: c fica na dimensao (nk-lb) x 1
                % -----------------------------------------------------------
                for ikprime = floor:nk
                    %%% Computar o consumo
                    c(ikprime)  = z.*(k^alpha) + (1-delta).*k - kgrid(ikprime);
                    
                    %%% Verificar consumos negativos
                    if c(ikprime) >= 0
                        u0 = (c(ikprime)^(1-mu)-1)/(1-mu);                        
                    else
                        u0 = -Inf;
                    end
                    
                    %%% Calcular a esperança (bem mais rápido)
                    Ev = 0;
                    for jz = 1:nz
                        Ev = Ev + P(iz,jz)*v(ikprime,jz);
                    end
                    
                    %%% Computar a função valor
                    H(ikprime) = u0 + beta.*Ev;
                    
                end
                
                %%% Computar a funçao valor (matriz)
                [Tv(ik,iz), idx_aux] = max(H(floor:nk));
                
                % Ajustamos a localização usando o lowerbound
                idx(ik,iz) = floor + idx_aux - 1;
                
                % Atualizamos o novo piso pra próxima iteração
                floor = idx(ik,iz);
            end
        end
        
    else
        
        for iz = 1:nz
            
            z = zgrid(iz);
            
            for ik = 1:nk
                
                % Pegando o capital ik e criando um vetor de consumo
                % -----------------------------------------------------------
                k = kgrid(ik);                  % por algum motivo aumenta a performance do cógido
                
                % Esse loop define o vetor de consumo de forma a garantir que
                % nenhuma entrada será negativa: a cada elemento, tiramos o
                % máximo entre o consumo e zero
                % -----------------------------------------------------------
                
                %%% Computar o consumo
                c(ik)  = z.*(k^alpha) + (1-delta).*k - kgrid(idx(ik,iz));
                
                %%% Verificar consumos negativos
                if c(ik) >= 0
                    u0 = (c(ik)^(1-mu)-1)/(1-mu);
                else
                    u0 = -Inf;
                end
                
                %%% Computar a esperança
                Ev = 0;
                for jz = 1:nz
                    Ev = Ev + P(iz,jz)*v(idx(ik,iz),jz);
                end
                
                %%% Computar a funçao valor (matriz)
                Tv(ik,iz) = u0 + beta.*Ev;
                
            end
        end
    end
    
    % Computa o erro e o número de iterações até agora
    error = max(max(abs(Tv - v)));
    iter  = iter + 1;
    
    % Atualiza o chute de V pelo novo TV encontrado
    v     = Tv;
    
    % Imprime a iteração e o erro
    fprintf('Error %4i %6.2e \n',[iter, error]);
    pause(0.001);
    clc;

end
time = toc;

    % Imprimir os resultados (tempo decorrido, nº de iterações e error)
    % -----------------------------------------------------------------
    last_iter = iter;

    if error<tol && last_iter<=max_iter

        fprintf('\n');
        fprintf('Method 2. Iteration of Value Function (Monotonicity w/ Accelerator) \n\n');
        fprintf('Time (in seconds)      = %3.03f \n',time);
        fprintf('Number of iterations   = %3.00f \n',last_iter);
        fprintf('Error                  = %3.10f \n',error);
        fprintf('\n');

    end
end