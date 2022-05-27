function [Tv, idx] = iteration_concavity_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

iter  = 0;
error = 1;

nk    = length(kgrid);
nz    = length(zgrid);

Tv    = zeros(nk,nz);
idx   = zeros(nk,nz);

% Iteração da Função valor - Concavidade c/ Acelerador
% ----------------------------------------------------
tic;
while  error > tol && iter <= max_iter
    
    if iter < 100 || rem(iter,10) == 0
       
        for iz = 1:nz
            
            z = zgrid(iz);
            
            for ik = 1:nk
                
                k = kgrid(ik);
                c = zeros(nk,1);
                H = zeros(nk,1);
                
                %%% Primeiro, computar para ikprime = 1
                ikprime = 1;
                
                %%% Computar o consumo em ikprime = 2
                c(ikprime) = z*(k^alpha) + (1-delta)*k - kgrid(ikprime);
                
                %%% Verificar consumos negativos
                if c(ikprime) >= 0
                    u0 = (c(ikprime)^(1-mu)-1)/(1-mu);
                else
                    u0 = -Inf;
                end
                
                %%% Calcular a esperança
                Ev = 0;
                for jz = 1:nz
                    Ev = Ev + P(iz,jz)*v(ikprime,jz);
                end
                
                %%% Computar a função valor
                H(ikprime) = u0 + beta.*Ev;
                
                %%% Agora, computar ikprime de 2 em diante
                for ikprime = 2:nk
                    
                    %%% Computar o consumo a partir de ikprime = 2
                    c(ikprime) = zgrid(iz)*(k^alpha) + (1-delta)*k - kgrid(ikprime);
                    
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
                    
                    % Calcular a convidade aqui! Veja que não há uso do max.
                    if H(ikprime) < H(ikprime-1)
                        break
                    end
                end
                
                %%% Mudança caso nunca tenha acionado o break (ou seja, o máximo é no último t mesmo)
                if  H(ikprime) <= H(ikprime-1)
                    
                    %%% Atualizar a função valor e a matriz de índices
                    Tv(ik,iz)  = H(ikprime-1); idx(ik,iz) = ikprime-1;
                    
                elseif H(ikprime) >= H(ikprime-1)
                    
                    %%% Atualizar a função valor e a matriz de índices
                    Tv(ik,iz)  = H(ikprime); idx(ik,iz) = ikprime;
                    
                end
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
    
    %%% Computar o erro
    error = max(max(abs(Tv - v)));
    iter = iter + 1;                    % atualizamos o contador de passos da iteração
    v = Tv;                             % atualizamos o chute de V pelo TV encontrado
    
    %%% Printar o número de iterações e o erro
    fprintf('Error %4i %6.2e \n', [iter, error]);
    pause(0.001);
    clc;

end
time = toc;

% Imprimir os resultados (tempo decorrido, nº de iterações e error)
% -----------------------------------------------------------------
last_iter = iter;

if error<tol && last_iter<=max_iter
    
    fprintf('\n');
    fprintf('Method 3. Iteration of Value Function (Concavity with Accelerator) \n\n');
    fprintf('Time (in seconds)      = %3.03f \n',time);
    fprintf('Number of iterations   = %3.00f \n',last_iter);
    fprintf('Error                  = %3.10f \n',error);
    fprintf('\n');
    
end
end