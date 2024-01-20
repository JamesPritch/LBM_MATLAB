%% Plotting of Benchmark 1a: advection-diffusion of a Gaussian hill
%% Finding analytical solution
% Load vars from simulation 
x = readmatrix("1a_x.txt");
y = readmatrix("1a_y.txt");
ux = readmatrix("1a_ux.txt");
uy = readmatrix("1a_uy.txt");
C_BGK = readmatrix("1a_c.txt");

% Set variables for analytical solution
x_0 = 200;
y_0 = 200;
t = 200;
D = 1.5;
C_0 = 1;
omega_0sqr = 10^2;
omega_Dsqr = 2*D*t;

% Find analytical solution for C
omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
xsqr = (x - x_0 - ux * t).^2;
ysqr = transpose((y - y_0 - uy * t).^2);
C_exact = omega_term * C_0 * exp( - (xsqr + ysqr) / (2 * (omega_0sqr+omega_Dsqr)));

% Error in C
Error = norm(C_BGK - C_exact, 'fro')/norm(C_exact, 'fro');


%% Plotting
% Plotting C against x
figure;
axes('FontSize',18, 'NextPlot', 'add');
plot(x, C_exact(200,:), 'k', x,C_BGK(200,:), 'k--', LineWidth=1)
% title("Concentration profile 1a", FontSize=20)
xlabel('${x/\Delta x}$','interpreter','latex', fontsize=30) 
ylabel('${C/C_0}$','interpreter','latex', fontsize=30)
xlim([100 300])
ylim([0 0.16])
leg = legend('Analytical', 'Numerical', 'interpreter','latex', fontsize = 22);
leg.ItemTokenSize = [25,25,25];
box on

% Plotting concentration contour plot 
figure;
axes('FontSize',18, 'NextPlot', 'add');
plot(1:2,1:2,'k');
plot(1:2,1:2,'k--');
contour(C_exact, [0.02 0.05 0.08 0.11 0.143], 'k', LineWidth=1)
[C,h] = contour(C_BGK, [0.02 0.05 0.08 0.11 0.145], 'k--', LineWidth=1);
clabel(C,h, 'FontSize',19, 'LabelSpacing',1000)
xlabel('${x/\Delta x}$','interpreter','latex', fontsize=30) 
ylabel('${y/\Delta y}$','interpreter','latex', fontsize=30)
axis equal
xlim([140 260])
ylim([140 260])
leg = legend({'Analytical', 'Numerical', '', ''}, 'interpreter','latex', ...
    fontsize = 22);
leg.ItemTokenSize = [25,25,25];
box on






