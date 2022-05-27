% discret.m
function v = discret(shocks, S, Pi)

    % Shocks: vetor de choques do processo contínuo; 
    % S: Espaço de Estados 
    % Pi: Matriz de transição
    
    % Criamos um vetor que registrará os estados do espaço discreto, com
    % base nos choques contínuos. 
    v = zeros(length(shocks),1)';
    
    % O índice de estado inicial é simplesmente a mediana, pois o grid é
    % centrado na média e tem # elementos ímpar.
    estado = find(S==median(S(:))); 
    
    % Para cada período dos 10,0000 simulados 
    for i = 1:length(shocks)
        % Registra a probabilidade de ir para o primeiro estado 
        cdf = Pi(estado, 1);  
        % Percorre cada um dos estados do grid discreto 
        for j = 1:length(S)
            % Se a CDF do choque é menor do que a probabilidade acumulada
            % até o estado anterior, registra transição para o estado
            % anterior
            if normcdf(shocks(i), 0, 0.007) < cdf
                % Atualizo o estado atual
                estado = j; 
         break
            else
                % Se a CDF do choque é menor do que a probabilidade acumulada até o
                % estado anterior, registra transição para o estado anterior 
                cdf = cdf + Pi(estado, j+1);
            end
        end
        v(i) = S(estado);
    end
end