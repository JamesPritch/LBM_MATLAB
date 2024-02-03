%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/model/T1
T = T(2:65,:); % This removes the upper and lower rows which weren't
               % part of the actual model

%% Setting Variables
% Scalars
nx = 234;
ny = 64;
scale = 1; % 18 is closest to mm in paper
a = 202; % Centre of curve at fingertip

% Vectors
x = 0:nx-1;
y = (0:ny-1) - (ny-1)/2;


%% Plotting T
% Temperature contour
figure;
axes('FontSize',11.5, 'NextPlot', 'add');
% Plotting temperature of finger
% surf(scale*x/(ny-1),scale*y/(ny-1),T);
s = pcolor(scale*x/(ny-1),scale*y/(ny-1),T);
s.LineWidth = 1e-10;
% Plotting finger outline
plot(zeros(ny),((0:ny-1) - (ny-1)/2)/(ny-1), 'k', LineWidth=1.9);
plot((0:a)/(ny-1),zeros(a+1) - 0.5, 'k', LineWidth=1.9);
plot((0:a)/(ny-1),zeros(a+1) + 0.5, 'k', LineWidth=1.9);
plot((a:nx)/(ny-1),-sqrt((1/2)^2 - (((a:nx) - a)/(ny-1)).^2), 'k', LineWidth=1.9);
plot((a:nx)/(ny-1),sqrt((1/2)^2 - (((a:nx) - a)/(ny-1)).^2), 'k', LineWidth=1.9);
% Making plot look nice
xlabel('${x/N_y}$','interpreter','latex', fontsize=19)
ylabel('${y/N_y}$','interpreter','latex', fontsize=19)
axis equal
xlim([-0.2*scale 3.85*scale])
ylim([-0.7*scale 0.7*scale])
colorbar
box on

% contour(scale*x/(ny-1),scale*y/(ny-1),outline, 'k')
% [C,h] = contour(scale*x/(ny-1),scale*y/(ny-1),T, 0.45:0.01:0.55, ...
%                 'k--', LineWidth=0.9);
% clabel(C,h, 'FontSize',12, 'LabelSpacing',502)



% Temperature cross section
% T_plot = (T(32,1:nx)-T_c)/(T_h-T_c);
% figure;
% axes('FontSize',11.5, 'NextPlot', 'add');
% plot(x/(ny-1),1-x/(2*(ny-1)),'k', ...
%      x/(ny-1),T_plot,'ok', LineWidth=0.9);
% xlabel('${x/N_x}$','interpreter','latex', fontsize=19)
% ylabel('${\frac{T-T_c}{T_h-T_c}}$','interpreter','latex', fontsize=23)
% axis equal
% xlim([0 2])
% ylim([-7e-4 1])
% leg = legend({'Expected', 'Numerical'}, 'interpreter','latex', ...
%     fontsize = 14);
% leg.ItemTokenSize = [12,12,12];
% box on











