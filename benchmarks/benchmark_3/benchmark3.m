%% Lattice Boltzmann method code: diffusion into steady flow
%% Setting Grid Independant Variables
% Simulation parameters - input
% Scalars
C_p = 1;
u_0 = 0.05;
nx = 1600;
ny = 160;
epsilon = 1e-7;
stop = false;

% Vectors
x = 1:nx;
y = 1:ny;

% D2Q9 velocity set parameters
ndir = 9;
cssq = 1/3;
cx = [0, 1, -1, 0, 0, 1, -1, -1, 1];
cy = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
cssqinv = 1/cssq;
omega = 1.38; % 1/tau;
omomega = 1 - omega;
D = 0.07407; % cssq*(tau - 0.5);

% Initialisation of the concentration and velocity field at time t = 0
C = zeros(ny, nx);
C(1,:) = C_p;

ux = zeros(ny, nx) + u_0;
uy = zeros(ny, nx);


%% Simulating using LBM
% Initialisation of the particle distribution function
geq = zeros(ny, nx, ndir);
for k = 1:ndir
    cdotu = cx(k)*ux + cy(k)*uy;
    udotu = ux.^2 + uy.^2;
    geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
end
g = geq;
gcol = zeros(ny, nx, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
while stop == false
    % Collision
    gcol = omomega*g + omega*geq;

    % Streaming - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for i = 1:ny
            for j = 1:nx
                xstreamed = mod(j + cx(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(i + cy(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(ystreamed, xstreamed, k) = gcol(i, j, k);
            end
        end
    end

    % Boundary conditions
    % LHS Boundary ie C = 0
    for i = 1:ny
        g(i,1,:) = 0;
    end
    % RHS Boundary ie partial C wrt x = 0
    for i = 1:ny
        g(i,nx,:) = g(i,nx-1,:);
    end
    % Top Boundary ie partial C wrt y = 0
    for i = 1:nx
        g(ny,i,:) = g(ny-1,i,:);
    end
    % Bottom Boundary ie Anti-BB
    opp = [0 2 1 4 3 7 8 5 6];
    for i = 1
        for j = 1:nx
            for k = 1:ndir
                g(i,j,opp(k)+1) = -gcol(i,j,k) + 2*w(k) * C_p;
            end
        end
    end

    % Macroscopic variables
    C1 = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
        + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);

    % Difference from past iteration
    diff = norm(C1-C)/norm(C);
    fprintf('Difference: %d, Time: %f \n', diff, toc);
    if diff < epsilon
        stop = true;
    end
    C = C1;

    % Equilibrium distribution function
    for k = 1:ndir
        cdotu = cx(k)*ux + cy(k)*uy;
        udotu = ux.^2 + uy.^2;
        geq(:, :, k) = w(k)*C.*(1 + cssqinv*cdotu + 0.5*cssqinv^2*cdotu.^2 - 0.5*cssqinv*udotu);
    end

end


%% Saving C to file:
save /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_3/C C









