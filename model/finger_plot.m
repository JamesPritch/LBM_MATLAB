%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/model/T1


%% Setting Variables
% Scalars
nx = 234;
ny = 64;
T_0 = 0.5; % T<1 seems to work
scale = 18;
T_c = T_0 - T_0/10; 
T_h = T_0 + T_0/10;
a = 202; % Centre of curve at fingertip

% Vectors
x = 0:nx-1;
y = (0:ny-1) - ny/2;
xs = zeros(ny,nx);
for i = 1:nx
    xs(:,i) = x(i);
end


%% Finger outline
outline = zeros(ny,nx);
outline( 1,1:a+1) = 1;
outline(ny,1:a+1) = 1;
r = zeros(ny,nx);
for i = 1:ny
    for j = 1:nx
        r(i,j) = sqrt((x(j)-a)^2+(y(i)+0.5)^2);
    end
end
[row,col] = find(abs(r)>31 & xs > a);
for i = 1:length(row)
    outline(row(i), col(i)) = 1;
end


%% Plotting T
% Temperature contour
figure;
axes('FontSize',11.5, 'NextPlot', 'add');
% plot(1:2,1:2,'k--');
contour(scale*x/(ny-1),scale*y/(ny-1),outline, 'k')
[C,h] = contour(scale*x/(ny-1),scale*y/(ny-1),T(1:ny,1:nx), 0.45:0.01:0.55, ...
                'k--', LineWidth=0.9);
clabel(C,h, 'FontSize',12, 'LabelSpacing',502)
xlabel('${x/N_x}$','interpreter','latex', fontsize=19)
ylabel('${y/N_y}$','interpreter','latex', fontsize=19)
axis equal
xlim([-0.2*scale 3.85*scale])
ylim([-0.7*scale 0.7*scale])
% leg = legend({'Numerical', ''}, 'interpreter','latex', ...
%     fontsize = 14);
% leg.ItemTokenSize = [25,25,25];
box on

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











