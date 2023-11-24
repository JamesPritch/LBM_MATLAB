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

    % Finding L_2 error
    Error(i) = norm(C_BGK{i} - C_exact{i}, 'fro')/norm(C_exact{i}, 'fro');

end


%% Plotting
x = 1:2000;
figure;
loglog([32 64 128 256 512], [Error(1) Error(2) Error(3) Error(4) Error(5)], ...
    'ko', x,4*x.^-1, 'k:', x,1000*x.^-2, 'k')
xlabel('Side Length of Square Domain (l.u.)','interpreter','latex', fontsize=26) 
ylabel('${L_2}$ error in ${C}$','interpreter','latex', fontsize=26)
xlim([10 2000])
ylim([10e-5 10e-0])
leg = legend('Simulations', 'First order Convergence', 'Second order Convergence', ...
             'interpreter','latex', fontsize = 16);
leg.ItemTokenSize = [12,12,12];






