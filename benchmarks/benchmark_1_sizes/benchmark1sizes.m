%% Lattice Boltzmann method code: advection-diffusion of a Gaussian hill
%% Setting Grid Independant Variables
% Simulation parameters - input
grids = 5;
niter_all = [3 12 48 192 768];
nx_all = [32 64 128 256 512];
ny_all = [32 64 128 256 512];
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
omega_0 = 10;
C_0 = 1;
x_0_all = [16 32 64 128 256];
y_0_all = [16 32 64 128 256];
ux_all = [0.2 0.1 0.05 0.025 0.0125];
uy_all = [0.2 0.1 0.05 0.025 0.0125];
C = cell(grids, 1);


%% Repeating LBM algorithm over grid sizes
for l=1:grids
    % Simulation parameters - input
    niter = niter_all(l);                                                               
    nx = nx_all(l);
    ny = nx_all(l);
    omega_0sqr = omega_0^2 * (nx_all(l)/nx_all(grids))^2;
    
    % Initialisation of the concentration field at time t = 0
    x = 1:nx;
    y = 1:ny;
    x_0 = zeros(nx, ny) + x_0_all(l);
    y_0 = zeros(nx, ny) + x_0_all(l);
    
    squaredterm = (x - x_0 + zeros(nx, ny)).^2 + transpose((y - y_0 + zeros(nx, ny)).^2);
    C{l} = C_0 * exp( - squaredterm / (2 * omega_0sqr));
    
    % Initialisation of the velocity field at time t = 0
    ux = zeros(nx, ny) + ux_all(l);
    uy = zeros(nx, ny) + ux_all(l);
    
    % Initialisation of the particle distribution function
    geq = zeros(nx, ny, ndir);
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C{l}.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end
    g = geq;
    gcol = zeros(nx, ny, ndir);
    
    % Simulation loop
    fprintf('Starting simulation %1.0f \n', l);
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
        C{l} = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
            + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    
        % Equilibrium distribution function
        for k = 1:ndir
            cdotu = cx(k)*ux + cy(k)*uy;
            udotu = ux.^2 + uy.^2;
            geq(:, :, k) = w(k)*C{l}.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
        end
    
        if mod(t, i^2) == 0
            fprintf('Iteration: %d, Time: %f \n', t, toc);
        end
    
    end

end


%% Saving C to file
save /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_1_sizes/C C



