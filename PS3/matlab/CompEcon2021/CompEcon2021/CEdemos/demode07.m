%% DEMODE07 Commercial Fisheries Model (V.L. Smith) 
%
%  Intertemporal price equilibrium in market for renewable resource.

%  Solve
%    s' = (1-s)*s - s*k/(alpha+beta*s*k) 
%    k' = delta*(alpha*s/(2(alpha+beta*s*k)^2}-phi)
%  where
%    s: Stock of fish
%    k: Industry size
%    x = [s;k]

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'Commercial Fisheries Model Phase Diagram';
x1label  = 'Fish Population';
x2label  = 'Industry Size';

% Velocity Function
alpha = 0.5;    % marginal cost parameter
beta  = 2.75;   % slope of demand curve
phi   = 0.05;   % fixed cost
delta = 10;     % industry entry rate
s = @(x) x(1,:);
k = @(x) x(2,:);
y = @(x) s(x)./(alpha+beta*s(x).*k(x));
p = @(x) alpha./(alpha+beta*s(x).*k(x));
f = @(x) [(1-s(x)).*s(x)-k(x).*y(x); delta*(0.5*p(x).*y(x)-phi)];

% Initial States
Ns = 7;
xinit = [nodeunif(Ns,0.0,1.0)';zeros(1,Ns)];
       
% Time Horizon
T = 30;


%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Time Discretization
N = 500;                % number of time nodes
       
% Solve for Different Initial Values
x  = zeros(N,2,Ns);
for i=1:Ns
  [x(:,:,i),t] = oderk4(f,xinit(:,i),T,N);
end

% Plot Solutions in Time Domain
for i=1:Ns
  figure
  plot(t,x(:,:,i))
  if i>1, ytickformat('%.1f'), end
  legend(x1label,x2label)
end


%% STEADY-STATES

% Steady State A
A = [0.1;0.6]; 
A = broyden(f,A);
disp('Steady State A')
disp(A)
disp('Eigenvalues')
disp(eig(fdjac(f,A)))

% Steady State B
B = [0.3;0.9]; 
B = broyden(f,B);
disp('Steady State B')
disp(B)
disp('Eigenvalues')
disp(eig(fdjac(f,B)))

% Steady State C
C = [0.5;0.8]; 
C = broyden(f,C);
disp('Steady State C')
disp(C)
disp('Eigenvalues')
disp(eig(fdjac(f,C)))


%% PHASE DIAGRAM

% Plotting Limits
x1lim = [0 1];    % Fish population plotting limits
x2lim = [0 1];    % Industry size plotting limits

% Compute Separatrix
xspx = odespx(f,B,T);

% Compute Nullclines
s = nodeunif(100,x1lim(1),x1lim(2));
snull = [s alpha*(1-s)./(1-beta*s.*(1-s))]; 
knull = [s (sqrt(alpha*s/(2*phi))-alpha)./(beta*s)];

% Plot Phase Diagram
odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,xspx,[A B C],snull,knull,1)
xtickformat('%.1f')
ytickformat('%.1f')

% Label Steady States
text(A(1)+0.03,A(2)-0.01,'A')
text(B(1)+0.01,B(2)-0.07,'B')
text(C(1)-0.04,C(2)-0.06,'C')


%% SAVE FIGURES
printfigures(mfilename)