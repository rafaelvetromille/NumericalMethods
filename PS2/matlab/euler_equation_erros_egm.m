function EEE = euler_equation_erros_egm(C1, g_aprox, kgrid, zgrid, P, alpha, beta, delta, mu)
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

nk    = length(kgrid);
nz    = length(zgrid);

EEE = zeros(nk,nz);

for iz = 1:nz
    for ik = 1:nk
        Ec = 0;
        for w = 1:nz
            
            aux1 = (interp1(kgrid, C1(:,w), g_aprox(ik,iz), 'linear'))^(-mu);
            aux2 = (zgrid(w)*alpha*(g_aprox(ik,iz)^(alpha-1)) + (1-delta));
            aux3 = aux1*aux2;
            
            Ec = Ec + beta*P(iz,w)*aux3;
        end
        EEE(ik,iz) = log10(abs(1 - (Ec^(-1/mu))/C1(ik,iz)));
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
title(['\textbf{Equa\c{c}{\~a}o dos Erros de Euler (EEE) para v{\''a}rios $$z$$''s}' newline '\textbf{(M{\''e}todo do Grid End{\''o}geno)}'],'interpreter','latex','fontsize',20);
ha = legend('$z=z_1$', '$z=z_2$', '$z=z_3$',...
    '$z=z_4$', '$z=z_5$', '$z=z_6$',...
    '$z=z_7$',...
    'Location','SouthEast');
set(ha,'Interpreter','latex');
axis([min(kgrid) max(kgrid) 1.05*min(min(EEE)) 0.95*max(max(EEE))])
end