%% Plotting of Benchmark 1a: advection-diffusion of a Gaussian hill
%% Finding analytical solution
% Load vars from simulation 
x = readmatrix("1b_x.txt");
y = readmatrix("1b_y.txt");
ux = readmatrix("1b_ux.txt");
uy = readmatrix("1b_uy.txt");
C_BGK = readmatrix("1b_c.txt");

% Set variables for analytical solution
x_0 = 200;
y_0 = 200;
t = 200;
D = 0.0043;
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
plot(x, C_exact(220,:), 'k', x,C_BGK(220,:), 'k--', LineWidth=1)
% title("Concentration profile 1b", FontSize=20)
xlabel('${x/\Delta x}$','interpreter','latex', fontsize=30) 
ylabel('${C/C_0}$','interpreter','latex', fontsize=30)
xlim([180 260])
leg = legend('Analytical', 'Numerical', 'interpreter','latex', fontsize = 22);
leg.ItemTokenSize = [25,25,25];
box on


% Plotting concentration contour plot 
figure;
axes('FontSize',18, 'NextPlot', 'add');
plot(1:2,1:2,'k');
plot(1:2,1:2,'k--');
contour(C_exact, [0.2 0.4 0.6 0.8], 'k', LineWidth=1)
[C,h] = contour(C_BGK, [0.2 0.4 0.6 0.8], 'k--', LineWidth=1);
clabel(C,h,'FontSize',19, 'LabelSpacing',1000)
% title("Concentration contour plot 1b", FontSize=20)
xlabel('${x/\Delta x}$','interpreter','latex', fontsize=30) 
ylabel('${y/\Delta y}$','interpreter','latex', fontsize=30)
axis equal
xlim([195 245])
ylim([195 245])
leg = legend({'Analytical', 'Numerical', '', ''}, 'interpreter','latex', ...
    fontsize = 22);
leg.ItemTokenSize = [25,25,25];
box on





