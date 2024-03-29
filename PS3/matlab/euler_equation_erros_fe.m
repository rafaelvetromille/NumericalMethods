function EEE = euler_equation_erros_fe(c_star, g_star, x, n_elements, parameters)

% Euler Equations Error (EEE)

% Derivada da fun��o de utilidade (caso geral)
if parameters.mu == 1
    u_prime  = @(c)(1./c);
else
    u_prime  = @(c)(c.^(-parameters.mu));
end

% Inv. da derivada da fun��o de utilidade (caso geral)
if parameters.mu == 1
    u_prime_inv  = @(c)(1./c);
else
    u_prime_inv  = @(c)(c.^(-1/parameters.mu));
end

% Computar os erros de Euler
nk    = length(parameters.kgrid);
nz    = length(parameters.zgrid);

EEE   = zeros(nk,nz);

for iz = 1:nz
    for ik = 1:nk
        
        aux1 = zeros(nz,1); aux2 = zeros(nz,1);
        
        for izz = 1:nz
            aux1(izz) = u_prime(c_fe(x(izz,:), g_star(ik,iz), n_elements, parameters));
            aux2(izz) = (parameters.zgrid(izz)*parameters.alpha*(g_star(ik,iz)^(parameters.alpha-1)) + (1-parameters.delta));
        end
        
        aux3  = aux1.*aux2;
        EEE(ik,iz) = log10(abs(1 - (u_prime_inv(parameters.beta*parameters.P(iz,:)*aux3))/c_star(ik,iz)));
        
    end
end

% Graph (2D) - Euler Equation Erros (EEE)
set(gcf,'DefaultLineLineWidth', 1.5);
hold on
for i = 1:nz
    plot(parameters.kgrid, EEE(:,i), 'DisplayName', sprintf('$z_{%4i} = %6.2f$', [i, parameters.zgrid(i)]));   
end
hl = legend('show', 'Location', 'southeast');
set(hl, 'Interpreter','latex')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('$\log10$ $\mid$ Euler Equation Error $\mid$','Interpreter','latex');
title(sprintf('\\textbf{Equa\\c{c}{\\~a}o dos Erros de Euler (EEE) para v{\\''a}rios $$z$$''s}\n\\textbf{(usando %4i pontos de coloca\\c{c}{\\~a}o + FE)}', n_elements), ...
    'Interpreter', 'latex', 'fontsize', 20)
legend('show')
axis([0.95*min(parameters.kgrid) 1.025*max(parameters.kgrid) 1.05*min(min(EEE)) 0.95*max(max(EEE))])
box on
hold off
end