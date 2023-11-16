%% Plotting of Benchamrk 1 with different grid sizes
%% Read C_BGK from file
load /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/C C
C_BGK = C;


%% Grid Independant Variables
% Scalars
grids = 5;
C_0 = 1;
D = 0.0043;
omega_0 = 10;

% Vectors
t = [3 12 48 192 768];
nx = [32 64 128 256 512];
ny = [32 64 128 256 512];
x_0 = [16 32 64 128 256];
y_0 = [16 32 64 128 256];
ux_all = [0.2 0.1 0.05 0.025 0.0125];
uy_all = [0.2 0.1 0.05 0.025 0.0125];

% Matrices
Error = zeros(1,grids);
C_exact = cell(grids, 1);


%% Loop to over grid sizes
for i=1:grids
    % Set variables for analytical solution
    x = 1:nx(i);
    y = 1:ny(i);
    omega_Dsqr = 2*D*t(i);
    omega_0sqr = omega_0^2 * (nx(i)/nx(grids))^2;
    ux = zeros(nx(i), ny(i)) + ux_all(i);
    uy = zeros(nx(i), ny(i)) + ux_all(i);

    % Find analytical solution for C
    omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
    xsqr = (x - x_0(i) - ux * t(i)).^2;
    ysqr = transpose((y - y_0(i) - uy * t(i)).^2);
    C_exact{i} = omega_term * C_0 * exp( - (xsqr + ysqr) / ...
        (2 * (omega_0sqr+omega_Dsqr)));

end


%% Plotting loop for different grid sizes
% Parameters to pre set
xpos_all = [17 33 66 133 266]; 
xlims = [[12 24 48  96 192]
         [22 44 88 176 352]];

% Plot each graph individually
for i=1:grids
    x = 1:nx(i);
    y = 1:ny(i);
    figure;
    axes('FontSize',16, 'NextPlot', 'add');
    plot(x, C_exact{i}(xpos_all(i),:), 'k', ...
         x, C_BGK{i}(xpos_all(i),:), 'k--')
    title(["Concentration profile Grid size ",num2str(16*2^i)], FontSize=20)
    xlabel('${x}$','interpreter','latex', fontsize=26) 
    ylabel('${C}$','interpreter','latex', fontsize=26)
    xlim([xlims(1,i) xlims(2,i)])
    legend('Analytical', 'BGK', fontsize = 18)
end



