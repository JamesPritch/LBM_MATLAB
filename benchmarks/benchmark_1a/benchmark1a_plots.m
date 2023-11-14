%% Plotting of Benchmark 1a: advection-diffusion of a Gaussian hill
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

% Test prints to ensure code running as expected
% print(xsqr)
% print(type(x))
% print(C_exact)
% print(type(C_exact))

% Plotting C against x
figure;
axes('FontSize',16, 'NextPlot', 'add');
plot(x, C_exact(200,:), 'k', x,C_BGK(200,:), 'k--')
title("Concentration profile 1a", FontSize=20)
xlabel('${x}$','interpreter','latex', fontsize=26) 
ylabel('${C}$','interpreter','latex', fontsize=26)
xlim([100 300])
legend('Analytical', 'BGK', fontsize = 18)


% Plotting concentration contour plot 
figure;
axes('FontSize',16, 'NextPlot', 'add');
contour(C_exact, [0.02 0.05 0.08 0.11 0.143], 'k', LineWidth=0.7)
contour(C_BGK, [0.02 0.05 0.08 0.11 0.145], 'k--', 'ShowText', 'on', ...
    LineWidth=0.7, LabelSpacing=1100)
title("Concentration contour plot 1a", FontSize=20)
xlabel('${x}$','interpreter','latex', fontsize=26) 
ylabel('${y}$','interpreter','latex', fontsize=26)
axis equal
xlim([140 260])
ylim([140 260])
legend('Analytical', 'BGK', fontsize = 18)






