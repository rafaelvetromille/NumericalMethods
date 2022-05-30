% Métodos Numéricos - EPGE/FGV 2022
% Professor: Cezar Santos
% Problem Set 3 - Aluno: Rafael Vetromille

% Início formal do documento
clearvars;
close all;
clc;

% Calibração
beta     = 0.987;                        % Fator de desconto; 
mu       = 2;                            % Coeficiente de aversão relativa ao risco;
alpha    = 1/3;                          % Share do capital na função de produção; 
delta    = 0.012;                        % Taxa de depreciação;        
const    = 0;                            % Intercepto do processo AR(1);
rho      = 0.95;                         % Parâmetro de persistência dos choques de renda;
sigma    = 0.007;                        % Variância do erro do processo de renda.

nz       = 7;                            % Número de estados exógenos para a versão discreta do processo de renda AR(1); 
m        = 3;                            % Número de desvios em relação à média do processo AR para definir a amplitude do espaço de estados. 

% Capital de estado estacionário
kss  = (1/alpha*(1/beta + delta - 1))^(1/(alpha-1));        % Capital de estado estacionário;
css  = kss^alpha + (1-delta)*kss - kss;                     % Consumo de estado estacionário;            

fprintf('Com a calibração proposta, o capital e o consumo de estado estacionário são %3.2f. e %3.2f, respectivamente. \n\n', [kss, css]);

% Método de Tauchen (já feito em listas anteriores)
[S,P] = mytauchen(const, m, rho, sigma, nz);

%% Exercício 1

% Para resolver esta lista, utilizaremos métodos de projeção. Para este
% item, resolva o problema usando um método de projeção global. Em
% particular, utilize polinômios de Chebyshev e o método da colocação
% (collocation points) para resolver o problema. Forneça evidência sobre
% sua solução: figuras da função valor e/ou função política, tempo de
% execução, Euler errors, etc.

% Construímos, primeiramente, o grid para o capital e o grid para a
% produtividade. Para o grid de capital, use 500 pontos linearmente espaçados
% no intervalo [0.75*kss, 1.25*kss]. 
nz       = 7; 
zmin     = -m*sqrt(sigma^2/(1 - rho^2));
zmax     = +m*sqrt(sigma^2/(1 - rho^2));

nk       = 500; 
kmax     = 1.25*kss;
kmin     = 0.75*kss;

kgrid    = linspace(kmin, kmax, nk)';
zgrid    = exp(linspace(zmin, zmax, nz))';

% Estrutura de parâmetros para passar nas funções
parameters.beta  = beta;
parameters.alpha = alpha;
parameters.delta = delta;
parameters.mu    = mu;
parameters.P     = P;
parameters.zgrid = zgrid;
parameters.kgrid = kgrid;

% Criamos uma função que cria os Polinômios de Chebychev e
% vamos testar essa função plotando os cinco primeiros polinômios:
% syms x
% fplot(chebyshev(0:4,x))
% axis([-1 1 -1 1])
% grid on
% ylabel('$T_n(x)$', 'Interpreter', 'latex')
% legend('$T_0(x)$','$T_1(x)$','$T_2(x)$','$T_3(x)$','$T_4(x)$','Location','northwest', 'Interpreter', 'latex')
% title('Chebyshev polynomials')

%%

% Turning off dialogs
options = optimset('Display', 'off');   
                   
% Otimização (atualizando o guess via 'for')
tic;

d = 5;      % Ordem do polinômico de Chebychev (-1)

for id = 1:d 
    if id == 1
        gamma0      = ones(nz, id+1);
        fun         = @(x) build_system(x, id, parameters);
        gamma_star  = fsolve(fun, gamma0, options);
        gamma0      = gamma_star;       % improve perfomance!!!
    else
        new_gamma0 = gamma0;
        for iz = 1:nz
            new_gamma0(iz,id+1) = 0;
        end
        fun        = @(x) build_system(x, id, parameters); 
        gamma_star = fsolve(fun, new_gamma0, options);
     end
end
timer(1) = toc;

%%

% Função consumo (c = c(k,z))
c_star = zeros(nk, nz);

% Preenchemos a função consumo usando o peso ótimo
for iz = 1:nz
    for ik = 1:nk
        c_star(ik,iz) = c_hat(gamma_star(iz,:), kgrid(ik), d, parameters);
    end
end

% Plotar a função política do consumo
figure(1)
plot_consumption_policy_function(c_star, kgrid, zgrid)

%%

% Função política do capital (k' = g(k,z))
g_star = zeros(nk, nz);

for iz = 1:nz
    for ik = 1:nk
        g_star(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c_star(ik,iz);
    end
end

% Plotar a função política do capital
figure(2)
plot_capital_policy_function(g_star, kgrid, zgrid)

%% Recuperando a função valor;
tic; 

% Recuperando a matriz de índices (está certo?)
idx = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        dif = abs(g_star(ik,iz) - kgrid);
        minimo = min(dif);
        idx(ik,iz) = find(minimo == dif);
    end
end

% Parâmetros numéricos
iter     = 0;
error    = 1;
v        = zeros(nk,nz);            % Guess inicial;
Tv       = zeros(nk,nz);            % Matriz do operador de Bellman;
tol      = 10e-6;
max_iter = 10000; 

while (error > tol && iter <= max_iter)
    for iz = 1:nz
        for ik = 1:nk
            
            %%% Pegando o capital ik e criando um vetor de consumo
            k = kgrid(ik);                  % por algum motivo aumenta a performance do cógido
    
            %%% Computar o consumo
            c         = zeros(nk,nz);
            c(ik,iz)  = zgrid(iz).*(k^alpha) + (1-delta).*k - kgrid(idx(ik,iz));
            
            %%% Computar a utilidade (não há consumo negativo)
            u0 = utility(c(ik,iz), mu);
            
            %%% Computar a esperança
            Ec = 0;
            for jz = 1:nz
                Ec = Ec + P(iz,jz)*v(idx(ik,iz),jz);
                %Ev = Ev + P(iz,jz)*interp1(kgrid, v(:,jz), g_star(ik,iz), 'linear');   % interpolando (mais correto?)
            end
            
            %%% Computar a funçao valor (matriz)
            Tv(ik,iz) = u0 + beta.*Ec;  
        end
    end
    
    % Calcular o erro e o número de iterações
    error = max(max(abs(Tv - v)));
    iter  = iter + 1;                 % atualiza a  iteração
    v     = Tv;                       % atualiza o chute de v pelo novo Tv encontrado
    
    % Imprime a iteração e o erro
    fprintf('Error %4i %6.2e \n', [iter, error]);
end
timer(2) = toc;

% Plotar a função valor
plot_value_function(v, kgrid, zgrid)

%% Erros de Euler (EEE)
EEE = euler_equation_erros(c_star, g_star, gamma_star, d, parameters);

%% Exercício 2

tic
nl  = 7;   % (L-1) intervalos
fun = @(x) build_system_fe(x, nl, parameters);

% Guess inicial
x0 = zeros(nz,nl);
for iz = 1:nz
    for il = 1:nl
        x0(iz,il) = il;
    end
end

% Otimização (atualizando o guess via 'for')
options = optimset('Display','off');     % Turning off dialogs
x_star  = fsolve(fun, x0, options);

timer = toc;

%%

% Função consumo (c = c(k,z))
c_new = zeros(nk,nz);

for i = 1:nz
    for j = 1:nk
        c_new(j,i) = c_fe(x_star(i,:), kgrid(j), nl, parameters);
    end
end

% Plotar a função política do consumo
figure(1)
plot_consumption_policy_function(c_new, kgrid, zgrid)

%%

% Função política do capital (k' = g(k,z))
g_new = zeros(nk, nz);

for iz = 1:nz
    for ik = 1:nk
        g_new(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c_new(ik,iz);
    end
end

% Plotar a função política do capital
figure(2)
plot_capital_policy_function(g_new, kgrid, zgrid)

%% Recuperando a função valor;
tic; 

% Recuperando a matriz de índices (está certo?)
idx = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        dif = abs(g_new(ik,iz) - kgrid);
        minimo = min(dif);
        idx(ik,iz) = find(minimo == dif);
    end
end

% Parâmetros numéricos
iter     = 0;
error    = 1;
v        = zeros(nk,nz);            % Guess inicial;
Tv       = zeros(nk,nz);            % Matriz do operador de Bellman;
tol      = 10e-6;
max_iter = 10000; 

while (error > tol && iter <= max_iter)
    for iz = 1:nz
        for ik = 1:nk
            
            %%% Pegando o capital ik e criando um vetor de consumo
            k = kgrid(ik);                  % por algum motivo aumenta a performance do cógido
    
            %%% Computar o consumo
            c         = zeros(nk,nz);
            c(ik,iz)  = zgrid(iz).*(k^alpha) + (1-delta).*k - kgrid(idx(ik,iz));
            
            %%% Computar a utilidade
            u0 = utility(c(ik,iz), mu);
            
            %%% Computar a esperança
            Ec = 0;
            for jz = 1:nz
                Ec = Ec + P(iz,jz)*v(idx(ik,iz),jz);
                %Ev = Ev + P(iz,jz)*interp1(kgrid, v(:,jz), g_star(ik,iz), 'linear');   % interpolando (mais correto?)
            end
            
            %%% Computar a funçao valor (matriz)
            Tv(ik,iz) = u0 + beta.*Ec;  
        end
    end
    
    % Calcular o erro e o número de iterações
    error = max(max(abs(Tv - v)));
    iter  = iter + 1;                 % atualiza a  iteração
    v     = Tv;                       % atualiza o chute de v pelo novo Tv encontrado
    
    % Imprime a iteração e o erro
    fprintf('Error %4i %6.2e \n', [iter, error]);
end
time = toc;

% Plotar a função valor
plot_value_function(v, kgrid, zgrid)

%% Erros de Euler (EEE)
EEE = euler_equation_erros_fe(c_new, g_new, x_star, nl, parameters);

