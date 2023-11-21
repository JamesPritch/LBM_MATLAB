%% Lattice Boltzmann method code: advection-diffusion of Cylinder without flow
%% Setting Grid Independant Variables
% Simulation parameters - input
% Scalars
grids = 2;
C_c = 1;
tau = 0.516;
a = 40;
nx = 129;
ny = 129;

% Vectors
niter = [16830 33692];
x = (1:nx) - 64.5;
y = (1:ny) - 64.5;
r = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        r(i,j) = sqrt(x(i)^2+y(j)^2);
    end
end

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
mask = zeros(nx, ny);
C = zeros(nx, ny, grids);
[row,col] = find(abs(r-a)<0.7);
for i = 1:length(row)
    mask(row(i), col(i)) = 1;
end
C(:,:,1) = mask(:,:);
C(:,:,2) = mask(:,:);

ux = zeros(nx, ny);
uy = zeros(nx, ny);

%% Simulating using LBM
% Initialisation of the particle distribution function
for l = 1:grids
    geq = zeros(nx, ny, ndir);
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C(:,:,l).*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end
    g = geq;
    gcol = zeros(nx, ny, ndir);
    
    % Simulation loop
    fprintf('Starting simulation \n');
    tic
    for t = 1:niter(l)
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
        opp = [0 2 1 4 3 7 8 5 6];
        for i = 1:length(row)
            for k = 1:ndir
                g(col(i),row(i),opp(k)+1) = -gcol(col(i),row(i),k) + 2*w(k) * C_c;
            end
        end
    
        % Macroscopic variables
        C(:,:,l) = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
            + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    
        % Equilibrium distribution function
        for k = 1:ndir
            cdotu = cx(k)*ux + cy(k)*uy;
            udotu = ux.^2 + uy.^2;
            geq(:, :, k) = w(k)*C(:,:,l).*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
        end
    
        if mod(t, 150) == 0
            fprintf('Iteration: %d, Time: %f \n', t, toc);
        end
    
    end

end

%% Saving C to file:
save /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_2/C C


