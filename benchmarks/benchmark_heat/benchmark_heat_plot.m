%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_heat/T



%% Setting Variables
% Scalars
% C_p = 1;
% D = 0.07407;
% u_0 = 0.05;
% nx = 1600;
% ny = 160;
% 
% % Vectors
% x = (1:nx);
% y = (1:ny);
% C_exact = zeros(ny,nx);


%% Finding Analytical Solution for C
% for i = 1:ny
%     for j = 1:nx
%         C_exact(i,j) = C_p * (1-erf(y(i)/sqrt(4*D*x(j)/u_0))) ;
%     end
% end


%% Plotting C contours
% Plotting C_exact and C_BGK Concentration contours 
% figure;
% axes('FontSize',11.5, 'NextPlot', 'add');
% plot(1:2,1:2,'k');
% plot(1:2,1:2,'k--');
% contour(x/ny,y/ny,C_exact, [0.2 0.4 0.6 0.8], 'k', LineWidth=0.9)
% [C,h] = contour(x/ny,y/ny,C_BGK, [0.2 0.4 0.6 0.8], 'k--', LineWidth=0.9);
% clabel(C,h,'manual', 'FontSize',12, 'LabelSpacing',502, 'Rotation',0)
% xlabel('${x/N_y}$','interpreter','latex', fontsize=19)
% ylabel('${y/N_y}$','interpreter','latex', fontsize=19)
% xlim([0 10])
% ylim([0 1])
% leg = legend({'Analytical', 'Numerical', '', ''}, 'interpreter','latex', ...
%     fontsize = 14);
% leg.ItemTokenSize = [25,25,25];
% box on

% Plotting C_BGK
figure;
contour(T)






