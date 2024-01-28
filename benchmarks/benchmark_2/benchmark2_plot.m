%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_2/C_antiBB
C_BGK_antiBB = C;
% load /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_2/C_inamuro 
% C_BGK_inamuro = C;


%% Setting Variables
% Scalars
C_c = 1;
D = 0.0052;
tau = 0.516;
a = 40;
nx = 129;
ny = 129;

% Vectors
t = [16830 33692];
x = (1:nx) - 64.5;
y = (1:ny) - 64.5;
ux = zeros(nx, ny);
uy = zeros(nx, ny);
r = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        r(i,j) = sqrt(x(i)^2+y(j)^2);
    end
end
mu = [2.4048 5.5201 8.6537 11.7915 14.9309 18.0703 21.2097 24.3491 27.4885];
C_exact = zeros(nx,ny,2);


%% Finding Analytical Solution for C
for k=1:2
    sumterm = 0;
    for i=1:length(mu)
        sumterm = sumterm + 2/(mu(i)*besselj(1,mu(i))) * ...
                    exp(-mu(i)^2*D*t(k)/a^2)*besselj(0,mu(i)*r/a);
    end
    C_exact(:,:,k) = C_c*(1-sumterm);
end


%% Plotting C against r/a
% Plotting parameters
xplot = 0:1/40:1;
xplot(1) = 1e-17;
C_plot_exact = zeros(2,41);
C_plot_exact(1,:) = C_exact(64,65:105,1);
C_plot_exact(2,:) = C_exact(64,65:105,2);
C_BGK_plot_antiBB = zeros(2,41);
C_BGK_plot_antiBB(1,:) = (C_BGK_antiBB(64,65:105,1) + C_BGK_antiBB(65,65:105,1))/2;
C_BGK_plot_antiBB(2,:) = (C_BGK_antiBB(64,65:105,2) + C_BGK_antiBB(65,65:105,2))/2;
% C_BGK_plot_inamuro = zeros(2,41);
% C_BGK_plot_inamuro(1,:) = (C_BGK_inamuro(64,65:105,1) + C_BGK_inamuro(65,65:105,1))/2;
% C_BGK_plot_inamuro(2,:) = (C_BGK_inamuro(64,65:105,2) + C_BGK_inamuro(65,65:105,2))/2;

% Error in C
Error = norm(C_BGK_plot_antiBB - C_plot_exact, 'fro')/norm(C_plot_exact, 'fro');

% Plotting C_exact and C_BGK Concentration profiles
figure;
axes('FontSize',14, 'NextPlot', 'add');
plot(xplot, C_plot_exact(1,:), 'k', LineWidth=1)
scatter(xplot, C_BGK_plot_antiBB(1,:), 50, 'diamondk', 'filled')
%scatter(xplot, C_BGK_plot_antiBB(1,:), 50, 'squarek', 'filled')
plot(xplot, C_plot_exact(2,:), ':k', LineWidth=1)
scatter(xplot, C_BGK_plot_antiBB(2,:), 50, 'diamondk')
%scatter(xplot, C_BGK_plot_antiBB(2,:), 50, 'squarek')
% title("Concentration profile for Benchmark 2", FontSize=20)
xlabel('${r/a}$','interpreter','latex', fontsize=23) 
ylabel('${C/C_c}$','interpreter','latex', fontsize=23)
xlim([0 1])
ylim([0 1])
leg = legend('${\hat{t} = 0.0547}$ (Analytical)', '${\hat{t} = 0.0547}$ (Numerical)', ...
    '${\hat{t} = 0.1095}$ (Analytical)', '${\hat{t} = 0.1095}$ (Numerical)', ...
    'interpreter','latex', 'Location','northwest', fontsize = 17);
%    '${\hat{t} = 0.1095}$ (Inamuro)', ...
%    '${\hat{t} = 0.0547}$ (Inamuro)', ...
leg.ItemTokenSize = [12,12,12];
box on


% Plotting C_BGK
% figure;
% contour(C_BGK_inamuro(:,:,1))

% Plotting mask
% figure;
% heatmap(mask)



