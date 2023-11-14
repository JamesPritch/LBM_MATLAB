%% Plotting of Benchamrk 1 with different grid sizes
%% Read variables from file
x_1024 = readmatrix("1b_x1024.txt");
y_1024 = readmatrix("1b_y1024.txt");
ux_1024 = readmatrix("1b_ux1024.txt");
uy_1024 = readmatrix("1b_uy1024.txt");
C_BGK_1024 = readmatrix("1b_c1024.txt");

x_512 = readmatrix("1b_x512.txt");
y_512 = readmatrix("1b_y512.txt");
ux_512 = readmatrix("1b_ux512.txt");
uy_512 = readmatrix("1b_uy512.txt");
C_BGK_512 = readmatrix("1b_c512.txt");

x_256 = readmatrix("1b_x256.txt");
y_256 = readmatrix("1b_y256.txt");
ux_256 = readmatrix("1b_ux256.txt");
uy_256 = readmatrix("1b_uy256.txt");
C_BGK_256 = readmatrix("1b_c256.txt");

x_128 = readmatrix("1b_x128.txt");
y_128 = readmatrix("1b_y128.txt");
ux_128 = readmatrix("1b_ux128.txt");
uy_128 = readmatrix("1b_uy128.txt");
C_BGK_128 = readmatrix("1b_c128.txt");


%% Grid size 512
% Set variables for analytical solution
x_0_1024 = 200;
y_0_1024 = 200;
t_1024 = 200;
D = 0.0043;
C_0 = 1;
omega_0sqr = 10^2;
omega_Dsqr = 2*D*t_512;

% Find analytical solution for C
omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
xsqr = (x_1024 - x_0_1024 - ux_1024 * t_1024).^2;
ysqr = transpose((y_1024 - y_0_1024 - uy_1024 * t_1024).^2);
C_exact_1024 = omega_term * C_0 * exp( - (xsqr + ysqr) / (2 * (omega_0sqr+omega_Dsqr)));

% Finding L_2 error
Error_1024 = sqrt(sum((C_BGK_1024-C_exact_1024).^2)/sum((C_exact_1024).^2));


%% Grid size 512
% Set variables for analytical solution
x_0_512 = 200;
y_0_512 = 200;
t_512 = 200;
D = 0.0043;
C_0 = 1;
omega_0sqr = 10^2;
omega_Dsqr = 2*D*t_512;

% Find analytical solution for C
omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
xsqr = (x_512 - x_0_512 - ux_512 * t_512).^2;
ysqr = transpose((y_512 - y_0_512 - uy_512 * t_512).^2);
C_exact_512 = omega_term * C_0 * exp( - (xsqr + ysqr) / (2 * (omega_0sqr+omega_Dsqr)));

% Finding L_2 error
Error_512 = sqrt(sum((C_BGK_512-C_exact_512).^2)/sum((C_exact_512).^2));


%% Grid size 256
% Set variables for analytical solution
x_0_256 = 100;
y_0_256 = 100;
t_256 = 800;
D = 0.0043;
C_0 = 1;
omega_0sqr = 10^2;
omega_Dsqr = 2*D*t_256;

% Find analytical solution for C
omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
xsqr = (x_256 - x_0_256 - ux_256 * t_256).^2;
ysqr = transpose((y_256 - y_0_256 - uy_256 * t_256).^2);
C_exact_256 = omega_term * C_0 * exp( - (xsqr + ysqr) / (2 * (omega_0sqr+omega_Dsqr)));

% Finding L_2 error
Error_256 = sqrt(sum((C_BGK_256-C_exact_256).^2)/sum((C_exact_256).^2));


%% Grid size 128
% Set variables for analytical solution
x_0_128 = 100;
y_0_128 = 100;
t_128 = 800;
D = 0.0043;
C_0 = 1;
omega_0sqr = 10^2;
omega_Dsqr = 2*D*t_128;

% Find analytical solution for C
omega_term = omega_0sqr / (omega_0sqr + omega_Dsqr);
xsqr = (x_128 - x_0_128 - ux_128 * t_128).^2;
ysqr = transpose((y_128 - y_0_128 - uy_128 * t_128).^2);
C_exact_128 = omega_term * C_0 * exp( - (xsqr + ysqr) / (2 * (omega_0sqr+omega_Dsqr)));

% Finding L_2 error
Error_128 = sqrt(sum((C_BGK_128-C_exact_128).^2)/sum((C_exact_128).^2));


%% Plotting
figure;
loglog([128 256 512 1024], [Error_128 Error_256 Error_512 Error_1024])
xlim([0 1000])







