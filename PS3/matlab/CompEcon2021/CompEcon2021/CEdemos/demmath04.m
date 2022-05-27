%% DEMMATH04 Standard Copulas
%
% Draws contour plots and scatter diagrams for Clayton and Gumbel copulas.

% Preliminary tasks
deminit(mfilename)

tau = 0.7;
nc = 5;
zbounds = 3;
n = 800;
u = linspace(0.01,0.99,n);
[U1,U2] = meshgrid(u,u);
U = [U1(:) U2(:)];

%% Plot Contours

figure

subplot(1,2,1)
hold on
c = copulapdf('Clayton',U,2*tau/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zbounds zbounds])
ylim([-zbounds zbounds])
plothdash([],0)
plotvdash(0,[])
ylim([-zbounds zbounds])
axis square
xticks([])
yticks([])
title('Clayton')
xlabel('$z_1$')
ylabel('$z_2$')

subplot(1,2,2)
hold on
c = copulapdf('Gumbel',U,1/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zbounds zbounds])
ylim([-zbounds zbounds])
plothdash([],0)
plotvdash(0,[])
axis square
xticks([])
yticks([])
title('Gumbel')
xlabel('$z_1$')
ylabel('$z_2$')


%% Scatter Plot

figure
n = 100;

subplot(1,2,1)
hold on
u = copularnd('Clayton',2*tau/(1-tau),n);
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zbounds zbounds])
ylim([-zbounds zbounds])
plothdash([],0)
plotvdash(0,[])
axis square
xticks([])
yticks([])
xlabel('$\tilde z_1$')
ylabel('$\tilde z_2$')
title('Clayton')

subplot(1,2,2)
hold on
u = copularnd('Gumbel',1/(1-tau),n);
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zbounds zbounds])
ylim([-zbounds zbounds])
plothdash([],0)
plotvdash(0,[])
axis square
xticks([])
yticks([])
title('Gumbel')
xlabel('$\tilde z_1$')
ylabel('$\tilde z_2$')


%% SAVE FIGURES
printfigures(mfilename)