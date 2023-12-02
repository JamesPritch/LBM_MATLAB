%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_3/C
C_BGK = C;


%% Setting Variables
% Scalars
C_p = 1;
D = 0.07407;
u_0 = 0.05;
nx = 1600;
ny = 160;

% Vectors
x = (1:nx);
y = (1:ny);
C_exact = zeros(ny,nx);


%% Finding Analytical Solution for C
for i = 1:ny
    for j = 1:nx
        C_exact(i,j) = C_p * (1-erf(y(i)/sqrt(4*D*x(j)/u_0))) ;
    end
end


%% Plotting C contours
% Plotting C_exact and C_BGK Concentration contours 
figure;
axes('FontSize',16, 'NextPlot', 'add');
plot(1:2,1:2,'k');
plot(1:2,1:2,'k--');
contour(x/ny,y/ny,C_exact, [0.2 0.4 0.6 0.8], 'k', LineWidth=0.7)
[C,h] = contour(x/ny,y/ny,C_BGK, [0.2 0.4 0.6 0.8], 'k--', LineWidth=0.7);
clabel(C,h,'FontSize',18, 'LabelSpacing',2000)
xlabel('${x/N_y}$','interpreter','latex', fontsize=26) 
ylabel('${y/N_y}$','interpreter','latex', fontsize=26)
xlim([0 10])
ylim([0 1])
leg = legend({'Analytical', 'BGK', '', ''}, 'interpreter','latex', ...
    fontsize = 18);
leg.ItemTokenSize = [25,25,25];

% Plotting C_BGK
% figure;
% contour(C_BGK)






