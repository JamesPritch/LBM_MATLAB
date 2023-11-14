%% Lattice Boltzmann method code: advection-diffusion of a Gaussian hill
%% Grid size 1024
% Simulation parameters - input
niter = 100;                                                               
nx = 1024;
ny = 1024;
tau = 0.513;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1/tau;
omomega = 1 - omega;

% Initialisation of the concentration and velocity field at time t = 0
x = 1:nx;
y = 1:ny;
x_0 = zeros(nx, ny) + 400;
y_0 = zeros(nx, ny) + 400;
omega_0 = 10;
C_0 = 1;

squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
C = C_0 * exp( - squaredterm / (2 * omega_0^2));

ux = zeros(nx, ny) + 0.2;
uy = zeros(nx, ny) + 0.2;

% Initialisation of the particle distribution function
geq = zeros(nx, ny, ndir);
for k = 1:ndir
    cdotu = cx(k)*ux + cy(k)*uy;
    udotu = ux.^2 + uy.^2;
    geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
end
g = geq;
gcol = zeros(nx, ny, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % Collision
    gcol = omomega*g + omega*geq;

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:ny
            for i = 1:nx
                xstreamed = mod(i + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(j + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(xstreamed, ystreamed, k) = gcol(i, j, k);
            end
        end
    end

    % Boundary conditions
    

    % Macroscopic variables
    C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
        + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end

    if mod(t, 10) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end

end

% Saving the objects:
writematrix(x,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_x1024.txt')
writematrix(y,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_y1024.txt')
writematrix(ux,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_ux1024.txt')
writematrix(uy,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_uy1024.txt')
writematrix(C,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_c1024.txt')


%% Grid size 512
% Simulation parameters - input
niter = 200;                                                               
nx = 512;
ny = 512;
tau = 0.513;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1/tau;
omomega = 1 - omega;

% Initialisation of the concentration and velocity field at time t = 0
x = 1:nx;
y = 1:ny;
x_0 = zeros(nx, ny) + 200;
y_0 = zeros(nx, ny) + 200;
omega_0 = 10;
C_0 = 1;

squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
C = C_0 * exp( - squaredterm / (2 * omega_0^2));

ux = zeros(nx, ny) + 0.1;
uy = zeros(nx, ny) + 0.1;

% Initialisation of the particle distribution function
geq = zeros(nx, ny, ndir);
for k = 1:ndir
    cdotu = cx(k)*ux + cy(k)*uy;
    udotu = ux.^2 + uy.^2;
    geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
end
g = geq;
gcol = zeros(nx, ny, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % Collision
    gcol = omomega*g + omega*geq;

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:ny
            for i = 1:nx
                xstreamed = mod(i + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(j + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(xstreamed, ystreamed, k) = gcol(i, j, k);
            end
        end
    end

    % Boundary conditions
    

    % Macroscopic variables
    C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
        + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end

    if mod(t, 10) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end

end

% Saving the objects:
writematrix(x,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_x512.txt')
writematrix(y,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_y512.txt')
writematrix(ux,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_ux512.txt')
writematrix(uy,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_uy512.txt')
writematrix(C,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_c512.txt')

%% Grid Size 256
% Simulation parameters - input
niter = 800;                                                               
nx = 256;
ny = 256;
tau = 0.513;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1/tau;
omomega = 1 - omega;

% Initialisation of the concentration and velocity field at time t = 0
x = 1:nx;
y = 1:ny;
x_0 = zeros(nx, ny) + 100;
y_0 = zeros(nx, ny) + 100;
omega_0 = 10;
C_0 = 1;

squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
C = C_0 * exp( - squaredterm / (2 * omega_0^2));

ux = zeros(nx, ny) + 0.05;
uy = zeros(nx, ny) + 0.05;

% Initialisation of the particle distribution function
geq = zeros(nx, ny, ndir);
for k = 1:ndir
    cdotu = cx(k)*ux + cy(k)*uy;
    udotu = ux.^2 + uy.^2;
    geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
end
g = geq;
gcol = zeros(nx, ny, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % Collision
    gcol = omomega*g + omega*geq;

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:ny
            for i = 1:nx
                xstreamed = mod(i + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(j + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(xstreamed, ystreamed, k) = gcol(i, j, k);
            end
        end
    end

    % Boundary conditions
    

    % Macroscopic variables
    C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
        + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end

    if mod(t, 10) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end

end

% Saving the objects:
writematrix(x,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_x256.txt')
writematrix(y,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_y256.txt')
writematrix(ux,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_ux256.txt')
writematrix(uy,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_uy256.txt')
writematrix(C,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_c256.txt')


%% Grid Size 128
% Simulation parameters - input
niter = 2400;                                                               
nx = 128;
ny = 128;
tau = 0.513;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1/tau;
omomega = 1 - omega;

% Initialisation of the concentration and velocity field at time t = 0
x = 1:nx;
y = 1:ny;
x_0 = zeros(nx, ny) + 50;
y_0 = zeros(nx, ny) + 50;
omega_0 = 10;
C_0 = 1;

squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
C = C_0 * exp( - squaredterm / (2 * omega_0^2));

ux = zeros(nx, ny) + 0.025;
uy = zeros(nx, ny) + 0.025;

% Initialisation of the particle distribution function
geq = zeros(nx, ny, ndir);
for k = 1:ndir
    cdotu = cx(k)*ux + cy(k)*uy;
    udotu = ux.^2 + uy.^2;
    geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
end
g = geq;
gcol = zeros(nx, ny, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % Collision
    gcol = omomega*g + omega*geq;

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for j = 1:ny
            for i = 1:nx
                xstreamed = mod(i + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(j + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(xstreamed, ystreamed, k) = gcol(i, j, k);
            end
        end
    end

    % Boundary conditions
    

    % Macroscopic variables
    C = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
        + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end

    if mod(t, 10) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end

end

% Saving the objects:
writematrix(x,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_x128.txt')
writematrix(y,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_y128.txt')
writematrix(ux,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_ux128.txt')
writematrix(uy,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_uy128.txt')
writematrix(C,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/1b_c128.txt')
