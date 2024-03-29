% M�todos Num�ricos - EPGE/FGV 2022
% Professor: Cezar Santos
% Problem Set 3 - Aluno: Rafael Vetromille

% In�cio formal do documento;
clearvars;
close all;
clc;

% Calibra��o;
beta     = 0.987;                        % Fator de desconto; 
mu       = 2;                            % Coeficiente de avers�o relativa ao risco;
alpha    = 1/3;                          % Share do capital na fun��o de produ��o; 
delta    = 0.012;                        % Taxa de deprecia��o;        
const    = 0;                            % Intercepto do processo AR(1);
rho      = 0.95;                         % Par�metro de persist�ncia dos choques de renda;
sigma    = 0.007;                        % Vari�ncia do erro do processo de renda.

nz       = 7;                            % N�mero de estados ex�genos para a vers�o discreta do processo de renda AR(1); 
m        = 3;                            % N�mero de desvios em rela��o � m�dia do processo AR para definir a amplitude do espa�o de estados. 

% Capital de estado estacion�rio;
kss  = (1/alpha*(1/beta + delta - 1))^(1/(alpha-1));        % Capital de estado estacion�rio;
css  = kss^alpha + (1-delta)*kss - kss;                     % Consumo de estado estacion�rio;            

fprintf('Com a calibra��o proposta, o capital e o consumo de estado estacion�rio s�o %3.2f e %3.2f, respectivamente. \n\n', [kss, css]);

% M�todo de Tauchen (j� feito em listas anteriores);
[S,P] = mytauchen(const, m, rho, sigma, nz);

%% Exerc�cio 1;

% Para resolver esta lista, utilizaremos m�todos de proje��o. Para este
% item, resolva o problema usando um m�todo de proje��o global. Em
% particular, utilize polin�mios de Chebyshev e o m�todo da coloca��o
% (collocation points) para resolver o problema. Forne�a evid�ncia sobre
% sua solu��o: figuras da fun��o valor e/ou fun��o pol�tica, tempo de
% execu��o, Euler errors, etc.

% Constru�mos, primeiramente, o grid para o capital e o grid para a
% produtividade. Para o grid de capital, use 500 pontos linearmente espa�ados
% no intervalo [0.75*kss, 1.25*kss]. 
nz       = 7; 
zmin     = -m*sqrt(sigma^2/(1 - rho^2));
zmax     = +m*sqrt(sigma^2/(1 - rho^2));

nk       = 500; 
kmax     = 1.25*kss;
kmin     = 0.75*kss;

kgrid    = linspace(kmin, kmax, nk)';
zgrid    = exp(linspace(zmin, zmax, nz))';

% Estrutura de par�metros para passar nas fun��es;
parameters.beta  = beta;
parameters.alpha = alpha;
parameters.delta = delta;
parameters.mu    = mu;
parameters.P     = P;
parameters.zgrid = zgrid;
parameters.kgrid = kgrid;

% Teste da fun��o que criei os Polin�mios de Chebychev;
syms x
fplot(chebyshev(0:4,x))
axis([-1 1 -1 1])
grid on
ylabel('$T_n(x)$', 'Interpreter', 'latex')
legend('$T_0(x)$','$T_1(x)$','$T_2(x)$','$T_3(x)$','$T_4(x)$','Location','northwest', 'Interpreter', 'latex')
title('Chebyshev polynomials')

%%

% Otimiza��o;
tic;

d = 5;      % Ordem do polin�mico de Chebychev (-1);

for id = 1:d 
    if id == 1
        gamma0      = ones(nz, id+1);
        fun         = @(x) build_system(x, id, parameters);
        gamma_star  = fsolve(fun, gamma0);
        gamma0      = gamma_star;       % Improve perfomance;
    else
        new_gamma0 = gamma0;
        for iz = 1:nz
            new_gamma0(iz,id+1) = 0;
        end
        fun        = @(x) build_system(x, id, parameters); 
        gamma_star = fsolve(fun, new_gamma0);
     end
end
timer(1) = toc;

%%

% Fun��o consumo (c = c(k,z));
c_star = zeros(nk, nz);
for iz = 1:nz
    for ik = 1:nk
        c_star(ik,iz) = c_hat(gamma_star(iz,:), kgrid(ik), d, parameters);
    end
end

% Plotar a fun��o pol�tica do consumo;
figure(1)
plot_consumption_policy_function(c_star, kgrid, zgrid)

%%

% Fun��o pol�tica do capital (k' = g(k,z));
g_star = zeros(nk, nz);
for iz = 1:nz
    for ik = 1:nk
        g_star(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c_star(ik,iz);
    end
end

% Plotar a fun��o pol�tica do capital;
figure(2)
plot_capital_policy_function(g_star, kgrid, zgrid)

%% Recuperando a fun��o valor;

% Recuperando a matriz de �ndices (est� certo?);
idx = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        dif = abs(g_star(ik,iz) - kgrid);
        minimo = min(dif);
        idx(ik,iz) = find(minimo == dif);
    end
end

tic; 
% Par�metros num�ricos;
iter     = 0;
error    = 1;
v        = zeros(nk,nz);            % Guess inicial;
Tv       = zeros(nk,nz);            % Matriz do operador de Bellman;
tol      = 10e-6;
max_iter = 10000; 

while (error > tol && iter <= max_iter)
    for iz = 1:nz
        for ik = 1:nk
            
            %%% Pegando o capital ik;
            k = kgrid(ik);                  % por algum motivo aumenta a performance do c�gido;
    
            %%% Computar o consumo;
            c         = zeros(nk,nz);
            c(ik,iz)  = zgrid(iz).*(k^alpha) + (1-delta).*k - kgrid(idx(ik,iz));
            
            %%% Computar a utilidade (n�o h� consumo negativo);
            u0 = utility(c(ik,iz), mu);
            
            %%% Computar a esperan�a;
            Ev = 0;
            for jz = 1:nz
                Ev = Ev + P(iz,jz)*v(idx(ik,iz),jz);
                %Ev = Ev + P(iz,jz)*interp1(kgrid, v(:,jz), g_star(ik,iz), 'linear');   % interpolando (mais correto?)
            end
            
            %%% Computar a fun�ao valor (matriz);
            Tv(ik,iz) = u0 + beta.*Ev;  
        end
    end
    
    % Calcular o erro e o n�mero de itera��es;
    error = max(max(abs(Tv-v)));
    iter  = iter + 1;                 % Atualiza a  itera��o;
    v     = Tv;                       % Atualiza o chute de v pelo novo Tv encontrado;
    
    % Imprime a itera��o e o erro;
    fprintf('Error %4i %6.2e \n', [iter, error]);
end
toc 

% Plotar a fun��o valor;
figure(3)
plot_value_function(v, kgrid, zgrid)

%% Erros de Euler (EEE);
EEE = euler_equation_erros(c_star, g_star, gamma_star, d, parameters);

mean(mean(EEE))

max(max(EEE))
min(min(EEE))
%% Exerc�cio 2.1

% Para este item, novamente utilize um m�todo de proje��o, mas, em vez
% de um m�todo espectral, use o m�todo dos elementos finitos. Divida o
% espa�o de estado em diversos elementos. Para resolver este problema,
% tente utilizar tanto o m�todo da coloca��o quanto Galerkin. De novo,
% evid�ncias!

tic;
n_elements = 11;   % (l-1) intervalos;
fun        = @(x) build_system_fe(x, n_elements, parameters);

% Guess inicial (sugerido pela monitora);
a0 = zeros(nz,n_elements);
for iz = 1:nz
    for il = 1:n_elements
        a0(iz,il) = il;
    end
end

% Otimiza��o;
options  = optimset('Display','off'); 
a_star   = fsolve(fun, a0, options);
timer(2) = toc;

%%

% Fun��o consumo (c = c(k,z));
c_new = zeros(nk,nz);
for i = 1:nz
    for j = 1:nk
        c_new(j,i) = c_fe(a_star(i,:), kgrid(j), n_elements, parameters);
    end
end

% Plotar a fun��o pol�tica do consumo;
figure(1)
plot_consumption_policy_function(c_new, kgrid, zgrid)

%%

% Fun��o pol�tica do capital (k' = g(k,z));
g_new = zeros(nk, nz);
for iz = 1:nz
    for ik = 1:nk
        g_new(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c_new(ik,iz);
    end
end

% Plotar a fun��o pol�tica do capital;
figure(2)
plot_capital_policy_function(g_new, kgrid, zgrid)

%% Recuperando a fun��o valor;

% Recuperando a matriz de �ndices (est� certo?);
idx = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        dif = abs(g_new(ik,iz) - kgrid);
        minimo = min(dif);
        idx(ik,iz) = find(minimo == dif);
    end
end

tic
% Par�metros num�ricos;
iter     = 0;
error    = 1;
v        = zeros(nk,nz);            % Guess inicial;
Tv       = zeros(nk,nz);            % Matriz do operador de Bellman;
tol      = 10e-6;
max_iter = 10000; 

while (error > tol && iter <= max_iter)
    for iz = 1:nz
        for ik = 1:nk
            
            %%% Pegando o capital ik;
            k = kgrid(ik);                  % Por algum motivo aumenta a performance do c�gido;
    
            %%% Computar o consumo;
            c         = zeros(nk,nz);
            c(ik,iz)  = zgrid(iz).*(k^alpha) + (1-delta).*k - kgrid(idx(ik,iz));
            
            %%% Computar a utilidade;
            u0 = utility(c(ik,iz), mu);
            
            %%% Computar a esperan�a;
            Ec = 0;
            for jz = 1:nz
                Ec = Ec + P(iz,jz)*v(idx(ik,iz),jz);
                %Ev = Ev + P(iz,jz)*interp1(kgrid, v(:,jz), g_star(ik,iz), 'linear');   % interpolando (mais correto?)
            end
            
            %%% Computar a fun�ao valor (matriz);
            Tv(ik,iz) = u0 + beta.*Ec;  
        end
    end
    
    % Calcular o erro e o n�mero de itera��es;
    error = max(max(abs(Tv-v)));
    iter  = iter + 1;                 % Atualiza a  itera��o;
    v     = Tv;                       % Atualiza o chute de v pelo novo Tv encontrado;
    
    % Imprime a itera��o e o erro;
    fprintf('Error %4i %6.2e \n', [iter, error]);
end
toc 

% Plotar a fun��o valor;
figure(3)
plot_value_function(v, kgrid, zgrid)

%% Erros de Euler (EEE);
EEE = euler_equation_erros_fe(c_new, g_new, a_star, n_elements, parameters);

max(max(EEE))

%% Exerc�cio 2.2

% Para este item, novamente utilize um m�todo de proje��o, mas, em vez
% de um m�todo espectral, use o m�todo dos elementos finitos. Divida o
% espa�o de estado em diversos elementos. Para resolver este problema,
% tente utilizar tanto o m�todo da coloca��o quanto Galerkin. De novo,
% evid�ncias!

tic;
n_elements  = 11;   % (l-1) intervalos;
n_int       = 10; 

fun = @(x) build_system_galerkin(x, n_elements, n_int, P, zgrid, kgrid, alpha, beta, delta, mu, nz);

% Guess inicial (sugerido pela monitora);
a0 = zeros(nz,n_elements);
for iz = 1:nz
    for il = 1:n_elements
        a0(iz,il) = il;
    end
end

% Otimiza��o;
options  = optimset('Display','off'); 
a_star1   = fsolve(fun, a0, options);
timer(2) = toc;

%%

% Fun��o consumo (c = c(k,z));
c_new1 = zeros(nk,nz);
for i = 1:nz
    for j = 1:nk
        c_new1(j,i) = c_fe(a_star1(i,:), kgrid(j), n_elements, parameters);
    end
end

% Plotar a fun��o pol�tica do consumo;
figure(1)
plot_consumption_policy_function(c_new1, kgrid, zgrid)

%%

% Fun��o pol�tica do capital (k' = g(k,z));
g_new1 = zeros(nk, nz);
for iz = 1:nz
    for ik = 1:nk
        g_new1(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c_new1(ik,iz);
    end
end

% Plotar a fun��o pol�tica do capital;
figure(2)
plot_capital_policy_function(g_new1, kgrid, zgrid)

%%

% Erros de Euler (Galerkin)
figure(3)
EEE = euler_equation_erros_galerkin(c_new1, g_new1, a_star1, n_elements, n_int, parameters);

max(max(EEE))

mean(mean(EEE))
