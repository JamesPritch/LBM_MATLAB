%% Plotting of Benchamrk 1 with different grid sizes
%% Read variables from file
load /Users/jpritch/Documents/MATLAB/benchmarks/steady_diffusion
T_numerical = T;

%% Setting Variables
% Scalars
nx = 60;
ny = 60;
T_c = 0.45; 
T_h = 0.55;
A = 4;

% Vectors
x = (0:nx-1)/(nx-1);
y = (0:ny-1)/(ny-1);


%% Analytical solution
T_exact = zeros(ny,nx);
for i = 1:ny
    T_exact(i,:) = (exp(A*y(i))-1)/(exp(A)-1);
end

T_numerical = (T_numerical - T_c)/(T_h-T_c);


%% Plotting T
% Initial setup
% figure;
% axes('FontSize',11.5, 'NextPlot', 'add');
% % Plotting finger outline
% plot(zeros(ny),((0:ny-1) - (ny-1)/2)/(ny-1), 'k', LineWidth=1);
% plot((0:a)/(ny-1),(zeros(a+1) - 0.5), 'k', LineWidth=1);
% plot((0:a)/(ny-1),zeros(a+1) + 0.5, 'k', LineWidth=1);
% th = linspace( pi/2, -pi/2, 100);
% R = 1/2;
% x_c = R*cos(th) + (a)/(ny-1);
% y_c = R*sin(th) + 0;
% plot(x_c,y_c, 'k', LineWidth=1);
% % Plotting bone outline
% plot((0:a)/(ny-1),zeros(a+1) - 0.25, 'k', LineWidth=0.7);
% plot((0:a)/(ny-1),zeros(a+1) + 0.25, 'k', LineWidth=0.7);
% R = 1/4;
% x_c = R*cos(th) + (a)/(ny-1);
% y_c = R*sin(th) + 0;
% plot(x_c,y_c, 'k', LineWidth=0.7);
% % Making plot look nice
% xlabel('${x/N_y}$','interpreter','latex', fontsize=19)
% ylabel('${y/N_y}$','interpreter','latex', fontsize=19)
% axis equal
% xlim([-0.2*scale 3.85*scale])
% ylim([-0.7*scale 0.7*scale])

% Concentration profile
figure;
axes('FontSize',18, 'NextPlot', 'add');
plot(x, T_exact(:,10), 'k', LineWidth=1)
plot(x, T_numerical(:,10), 'ko', LineWidth=1)
% title("Concentration profile 1a", FontSize=20)
xlabel('${y/\Delta y}$','interpreter','latex', fontsize=30) 
ylabel('${\frac{T-T_0}{T_h-T_c}}$','interpreter','latex', fontsize=30)
% xlim([0 1])
% ylim([0 1])
leg = legend('Analytical', 'Numerical', 'interpreter','latex', Location = 'northwest', fontsize = 22);
leg.ItemTokenSize = [25,25,25];
box on

% Temperature heatmap
% figure;
% axes('FontSize',18, 'NextPlot', 'add');
% % Plotting temperature of finger
% s = pcolor(x/(ny-1),y/(ny-1),T_numerical);
% set(s, 'EdgeColor', 'none');
% % Making plot look nice
% xlabel('${x/N_y}$','interpreter','latex', fontsize=25)
% ylabel('${y/N_y}$','interpreter','latex', fontsize=25)
% axis equal
% xlim([0 1])
% ylim([0 1])
% c = colorbar;
% c.Limits = [T_c T_h];
% box on

% Temperature contour
% figure;
% axes('FontSize',18, 'NextPlot', 'add');
% % Plotting temperature of finger
% contour(x/(ny-1),y/(ny-1),T_numerical, 0.45:0.005:0.55, ...
%                 LineWidth=0.9);
% % Making plot look nice
% xlabel('${x/N_y}$','interpreter','latex', fontsize=25)
% ylabel('${y/N_y}$','interpreter','latex', fontsize=25)
% axis equal
% xlim([0 1])
% ylim([0 1])
% c = colorbar;
% c.Limits = [T_c T_h];
% box on












