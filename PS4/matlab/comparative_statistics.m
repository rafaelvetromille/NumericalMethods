function z = comparative_statistics(mu, rho, sig)

% a) (Espaço de Estados Exógenos) Utilize o Método de Tauchen para discretizar o processo de renda. 
% Para tanto, use N=9 (número de pontos nogrid) e m=3 (m desvios em relação à média para cada lado). 
% Apresente a matriz de transição P resultante e o espaço de estados discreto S para as realizações de renda.

% Resposta:

% Calibração
const    = 0;                            % Intercept term of the AR(1) process
beta     = 0.96;                         % Fator de desconto;
phi      = 1;                            % Limite de endividamento; 

ny       = 9;                            % Número de estados exógenos para a versão discreta do processo de renda AR(1); 
m        = 3;                            % Número de desvios em relação à média do processo AR para definir a amplitude do espaço de estados. 

% Método de Tauchen

% FUNCTION: [s, Pi] = mytauchen(const,m,rho,sig,N)
%
% INPUTS:
%   const:  scalar, intercept term of the AR(1) process
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%   N:      scalar, number of grid points for the discretized process
%
% OUTPUTS:
%   S:      column vector of size Nx1, contains all possible states in ascending order
%   Pi:     matrix of size NxN, contains the transition proabilities. Rows are current state and columns future state

[S, P] = mytauchen(const, m, rho, sig, ny);

%% %% EXERCISE 1B

% b) Discretize o espaço de ativos numgridde 201pontos, estabelecendo amax=4. 
% Impondo r=0.04, resolva o problema das famílias para encontrar as funções políticas (consumo e ativos) e a função valor. 
% Use gráficos para ilustrar essas três funções.

% Resposta:

ymin       = -m*sqrt(sig^2/(1 - rho^2));
ymax       = +m*sqrt(sig^2/(1 - rho^2));

na         = 201; 
amin       = -phi;
amax       = 4; 

agrid      = nodeunif(na, amin, amax);  
ygrid      = exp(nodeunif(ny, ymin, ymax));

%%%%% numerical parameters

max_iter   = 1000;
tol        = 1e-7;
penalty    = 10^16;

%%%% state grid

s          = gridmake(agrid,ygrid); % ns-by-2 matrix where ns=na*ny
ns         = size(s,1);        

a          = s(:,1);
y          = s(:,2);

%%%%% parameters structure to pass to functions

parameters.beta    = beta;
parameters.mu      = mu;
parameters.agrid   = agrid;
parameters.a       = a;
parameters.y       = y;
parameters.P       = P;

% GIVEN R RECOVER OTHER VARIABLES OF INTEREST
r = 0.04;

%%%%% Reward Function

c  = zeros(ns,na);

for j=1:na
    
    c(:,j) = (1+r)*a + y - agrid(j);
  
end

% CHECK CONSUMPTION NON-NEGATIVE

violations = (c<0);

c = c.*(c>=0) + eps;  % soma um valor muito pequeno para não ficar zero

% UTILITY FUNCTION

if mu==1
    u  = @(c)log(c);
    u0 = log(c) - penalty*violations;
else    
    u  = @(c)(c.^(1-mu)-1)./(1-mu);
    u0 = (1/(1-mu))*(c.^(1-mu) - 1) - penalty*violations;
end

% INITIALIZE VALUE FUNCTION
Vguess = zeros(na,ny);

for iy = 1:ny
    Vguess(:,iy) = u(r.*agrid+ygrid(iy))./(1-beta);   % guess sugerido Benjamin Moll
end

v1       = log(0.5*y)/(1-beta);                       % initial guess 1
v2       = (1/(1-mu))*((r.*a+y).^(1-mu) - 1);         % initial guess 2
v3       = reshape(Vguess,[],1);                      % initial guess 3
v4       = nodeunif(1809, 0, 0);                      % initial guess 4

% ITERATION OF VALUE FUNCTION

v       = v3;  % choose a initial guess (all of them give us the same result)
iter    = 0;

for i=1:max_iter

RHS      = u0+beta*kron(P,ones(na,1))*reshape(v,na,9)';

[Tv,argmax] = max(RHS,[],2);

%%%%% policy that attains the maximum

g = a(argmax);

%%%%% check if converged

error = max(abs(Tv-v));

% fprintf('Error %4i %6.2e \n',[i, error]);

if norm(Tv-v,inf)<tol, break, end

%%%%% if not converged, update guess

v = Tv;

end

% Compute the consumption (long format)
c_long   = (1+r)*a+y-g;

% Value function: V(a,y)
v        = reshape(v,na,ny);

% Graph (2D)
% figure(1)
% set(gcf,'DefaultLineLineWidth', 1.5);
% plot(agrid,v(:,1),agrid,v(:,2),agrid,v(:,3),...
%      agrid,v(:,4),agrid,v(:,5),agrid,v(:,6),...
%      agrid,v(:,7),agrid,v(:,8),agrid,v(:,9))
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
% 
% xlabel('asset level, $a$','Interpreter','latex');
% ylabel('value function, $v(a,\bar{y})$','Interpreter','latex');
% ha = legend('$y=y_1$', '$y=y_2$', '$y=y_3$',...
%             '$y=y_4$', '$y=y_5$', '$y=y_6$',...
%             '$y=y_7$', '$y=y_8$', '$y=y_9$',...
%             'Location','SouthEast');
% set(ha,'Interpreter','latex');
% axis([amin amax min(min(v)) max(max(v))]) 

% Policy functions: c(a,y) and g(a, y)
c        = reshape(c_long,na,ny);

% figure(2)
% set(gcf,'DefaultLineLineWidth', 1.5);
% plot(agrid,c(:,1),agrid,c(:,2),agrid,c(:,3),...
%      agrid,c(:,4),agrid,c(:,5),agrid,c(:,6),...
%      agrid,c(:,7),agrid,c(:,8),agrid,c(:,9),...
%      agrid,ones(size(agrid)),'k--',ones(size(agrid)),agrid,'k--')
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
% 
% xlabel('asset level, $a$','Interpreter','latex');
% ylabel('consumption policy, $c(a,y)$','Interpreter','latex');
% ha = legend('$y=y_1$', '$y=y_2$', '$y=y_3$',...
%             '$y=y_4$', '$y=y_5$', '$y=y_6$',...
%             '$y=y_7$', '$y=y_8$', '$y=y_9$',...
%             'Location','SouthEast');
% set(ha,'Interpreter','latex');
% axis([amin amax min(min(c)) max(max(c))]) 

% Policy functions: a'=g(a,y)
g        = reshape(g,na,ny);

% figure(3)
% set(gcf,'DefaultLineLineWidth', 1.5);
% plot(agrid,g(:,1),agrid,g(:,2),agrid,g(:,3),...
%      agrid,g(:,4),agrid,g(:,5),agrid,g(:,6),...
%      agrid,g(:,7),agrid,g(:,8),agrid,g(:,9),...
%      agrid,ones(size(agrid)),'k--',ones(size(agrid)),agrid,'k--')
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
% hold off
% xlabel('asset level, $a$','Interpreter','latex');
% ylabel('asset policy, $a^{\prime}=g(a,y)$','Interpreter','latex');
% ha = legend('$y=y_1$', '$y=y_2$', '$y=y_3$',...
%             '$y=y_4$', '$y=y_5$', '$y=y_6$',...
%             '$y=y_7$', '$y=y_8$', '$y=y_9$',...
%             'Location','SouthEast');
% set(ha,'Interpreter','latex');
% axis([amin amax amin amax]) 

%% EXERCISE 1C - ~3min 

% c) (A distribuição invariante): Encontre a distribuição invariante conjunta sobre (a, s)
% usando um dos métodos discutidos em monitoria e comente. Calcule, represente em um
% gráfico e comente também a distribuição invariante sobre os ativos apenas (distribuição
% de riqueza).

%%%%% CONSTRUCT TRANSITION MATRIX FOR THE STATE S=(A,Y)

A = zeros(ns,na);
Q = zeros(ns,ns);

PP = kron(P,ones(na,1));

f = waitbar(0, 'Starting');

for s=1:ns
    
    A(s,:) = (agrid==g(s))';  %% puts a 1 if g(s)=a 
    
    Q(s,:) = kron(PP(s,:),A(s,:));
    
    waitbar(s/ns, f, sprintf('Progress: %d %%', floor(s/ns*100)));
    pause(0.1);
    
end

close(f)

%%%%% COMPUTE STATIONARY DISTRIBUTION

[eig_vectors,eig_values] = eig(Q'); 
[~,arg] = min(abs(diag(eig_values)-1)); 
unit_eig_vector = eig_vectors(:,arg); 

lambda = unit_eig_vector/sum(unit_eig_vector); 
psi = reshape(lambda, na, ny);

%%%%% graph - distribuição invariante sobre os ativos apenas (distribuição de riqueza).

% figure(3)
% set(gcf,'DefaultLineLineWidth', 1.5);
% plot(agrid,psi(:,1)+psi(:,2)+psi(:,3)+...
%            psi(:,4)+psi(:,5)+psi(:,6)+...
%            psi(:,7)+psi(:,8)+psi(:,9))
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
% xlabel('asset level, $a$','Interpreter','latex');
% ylabel('marginal asset distribution, $\sum_y \lambda(a,y)$','Interpreter','latex');
 
% figure(3)
% set(gcf,'DefaultLineLineWidth', 1.5);
% plot(agrid,sum(psi,2))
% set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
% set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')  
% xlabel('asset level, $a$','Interpreter','latex');
% ylabel('marginal asset distribution, $\sum_y \lambda(a,y)$','Interpreter','latex');

%% EXERCISE 1D

% d) Apure o excesso de demanda líquida pelo ativo calculado sob a taxa de juros proposta.

g = reshape(g, [], 1);
z = sum(lambda.*g);

end
