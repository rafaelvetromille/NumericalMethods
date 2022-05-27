function plot_value_function(variable, kgrid, zgrid)

% Necess�rio para plotar o gr�fico em 3D
k            = repelem(kgrid,1,7);
z            = repelem(zgrid,1,length(kgrid))';
nz           = length(zgrid);

% Gr�fico (3D) - Fun��o valor: V(k,z)
subplot(1,2,1);
set(gcf,'DefaultLineLineWidth', 1.5);
surf(k, z, variable, 'FaceAlpha', 0.8, 'EdgeColor', 'none')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex')
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('realiza\c{c}{\~o}es de produtividade, $z$','Interpreter','latex');
zlabel('fun\c{c}{\~a}o valor, $v(k,z)$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o valor $$v(k,z)$$}','interpreter','latex','fontsize',20);
axis([min(kgrid) max(kgrid) min(zgrid) max(zgrid) min(min(variable)) max(max(variable))]) 
box on

% Gr�fico (2D) - Fun��o valor: V(k,zbar)
subplot(1,2,2);
set(gcf,'DefaultLineLineWidth', 1.5);
hold on
for i = 1:nz
    plot(kgrid, variable(:,i), 'DisplayName', sprintf('$z_{%4i} = %6.2f$', [i, zgrid(i)]));   
end
legend('show', 'Location', 'southeast')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('fun\c{c}{\~a}o valor, $v(k,\bar{z})$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o valor $$v(k,z)$$ para v{\''a}rios $$z$$''s}','interpreter','latex','fontsize',20);
axis([min(kgrid) max(kgrid) min(min(variable)) max(max(variable))]) 
box on 

end