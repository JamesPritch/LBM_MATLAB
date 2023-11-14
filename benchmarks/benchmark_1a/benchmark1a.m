%% Lattice Boltzmann method code: advection-diffusion of a Gaussian hill
% Simulation parameters - input
niter = 200;                                                               
nx = 512;
ny = 512;
tau = 5;

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
D = cssq*(tau - 0.5);

% Initialisation of the concentration and velocity field at time t = 0
x = 1:nx;
y = 1:ny;
x_0 = zeros(nx, ny) + 200;
y_0 = zeros(nx, ny) + 200;
omega_0 = 10;
C_0 = 1;

squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
C = C_0 * exp( - squaredterm / (2 * omega_0^2));

ux = zeros(nx, ny);
uy = zeros(nx, ny);

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
writematrix(x,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1a/1a_x.txt')
writematrix(y,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1a/1a_y.txt')
writematrix(ux,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1a/1a_ux.txt')
writematrix(uy,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1a/1a_uy.txt')
writematrix(C,'/Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1a/1a_c.txt')



