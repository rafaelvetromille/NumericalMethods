function EEE = euler_equation_erros(idx, c, kgrid, zgrid, P, alpha, beta, delta, mu)

% Euler Equations Error (EEE)

% Derivada da função de utilidade (caso geral)
% --------------------------------------------
if mu == 1
    u_prime  = @(c)(1./c);
else
    u_prime  = @(c)(c.^(-mu));
end

% Inv. da derivada da função de utilidade (caso geral)
% ----------------------------------------------------
if mu == 1
    u_prime_inv  = @(c)(1./c);
else
    u_prime_inv  = @(c)(c.^(-1/mu));
end

% Computar os erros de Euler
% --------------------------

nk    = length(kgrid);
nz    = length(zgrid);

EEE = zeros(nk,nz);

for iz = 1:nz
    for ik = 1:nk
        
        index           = idx(idx(ik,iz),:)';
        
        aux1            = u_prime(zgrid*(kgrid(idx(ik,iz))^alpha) + (1-delta)*kgrid(idx(ik,iz)) - kgrid(index));
        aux2            = zgrid.*alpha.*kgrid(idx(ik,iz))^(alpha-1) + (1-delta);
        
        aux3            = aux1.*aux2;
        
        EEE(ik,iz)      = log10(abs(1-u_prime_inv(P(iz,:)*(beta.*aux3))/c(ik,iz)));
        
    end
end

% Graph (2D) - Euler Equation Erros (EEE) - Força Bruta
% -----------------------------------------------------
set(gcf,'DefaultLineLineWidth', 1.5);
plot(kgrid,EEE)
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')
hold off
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('$\log10$ $\mid$ Euler Equation Error $\mid$','Interpreter','latex');
title(sprintf('\\textbf{Equa\\c{c}{\\~a}o dos Erros de Euler (EEE) para v{\\''a}rios $$z$$''s}\n\\textbf{(grid do capital com %4i pontos)}', nk), ...
    'Interpreter','latex','fontsize',20)
ha = legend('$z=z_1$', '$z=z_2$', '$z=z_3$',...
    '$z=z_4$', '$z=z_5$', '$z=z_6$',...
    '$z=z_7$',...
    'Location','SouthEast');
set(ha,'Interpreter','latex');
axis([min(kgrid) max(kgrid) 1.05*min(min(EEE)) 0.95*max(max(EEE))])

end