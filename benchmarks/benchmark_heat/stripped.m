%% Lattice Boltzmann method code: convection in square cavity
%% Setting Grid Independant Variables
% Simulation parameters - input
% Scalars
nx = 64;
ny = 64;
niter = 1000;
T_0 = 0.5; % This is average between 0ºC and 20ºC ie 10ºC. likes <1
T_c = 0.45; % Set to non-dim and find proper non-dim later
T_h = 0.55; % doesn't seem to like any bigger than this but smaller is fine
rho_0 = 1.0000005;
beta = 1e-7; % Check hand calculation of this 100. likes <1e-2
grav = 1; % 9.8 in dim, idk for non-dim. likes <10

% Vectors
x = 0:nx-1;
y = 0:ny-1;

% Matrices
rho = zeros(ny,nx) + rho_0; % Density

% D2Q9 velocity set parameters
ndir = 9;
c_s = sqrt(1/3);
c = 1; % Assuming this is 1
zeta_x = [0, 1, -1, 0, 0, 1, -1, -1, 1];
zeta_y = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
tau_c = 1.24;
tau_v = 1; % I don't know this, this is a guess
v = c^2 * (tau_v - 0.5);
chi = (2/3) * c^2 * (tau_c - 0.5);

% Initialisation of the temperature and velocity field at time t = 0
T = zeros(ny, nx) + T_0;
T(:, 1) = T_c;
T(:,ny) = T_h;

ux = zeros(ny, nx);
uy = zeros(ny, nx);

% Force fields
G = zeros(ny,nx,ndir);
F = zeros(ny,nx,ndir);


%% Simulating using LBM
% Initialisation of the particle distribution function
% Fluid density f
feq = zeros(ny, nx, ndir);
for k = 1:ndir
    feq(:, :, k) = w(k)*rho;
end
f = feq;
fcol = zeros(ny, nx, ndir);

% Fluid temperature g
geq = zeros(ny, nx, ndir);
geq(:, :, 1) = 0;
for k = 2:5
    geq(:, :, k) = rho.* T / 9 * 1.5;
end
for k = 6:9
    geq(:, :, k) = rho.* T / 36 * 3;
end
g = geq;
gcol = zeros(ny, nx, ndir);

% Simulation loop
fprintf('Starting simulation \n');
tic
for t = 1:niter
    % LBM for fluid density f first
    % Collision
    fcol = f - f * 1/tau_v + feq * 1/tau_v + F;

    % Streaming for f - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for i = 1:ny
            for j = 1:nx
                xstreamed = mod(j + zeta_x(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(i + zeta_y(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                f(ystreamed, xstreamed, k) = fcol(i, j, k);
            end
        end
    end

    % Boundary conditions for f
    % LHS Boundary ie u = 0
    % for i = 1:ny
    %     for k = [2 6 9]
    %         f(i,1,k) = 0;
    %     end
    % end
    % RHS Boundary ie u = 0
    % for i = 1:ny
    %     for k = [4 7 8]
    %         f(i,nx,:) = 0;
    %     end
    % end
    % Top Boundary ie Anti-BB
    % opp = [0 2 1 4 3 7 8 5 6];
    % for i = ny
    %     for j = 1:nx
    %         for k = 1:ndir
    %             f(i,j,opp(k)+1) = -fcol(i,j,k) + 2*w(k) * C_p;
    %         end
    %     end
    % end
    % Bottom Boundary ie Anti-BB
    % for i = 1
    %     for j = 1:nx
    %         for k = 1:ndir
    %             f(i,j,opp(k)+1) = -fcol(i,j,k) + 2*w(k) * C_p;
    %         end
    %     end
    % end

    % Macroscopic variables from f
    rho = f(:, :, 1) + f(:, :, 2) + f(:, :, 3) + f(:, :, 4) + f(:, :, 5)...
        + f(:, :, 6) + f(:, :, 7) + f(:, :, 8) + f(:, :, 9);
    rho_x_ux = f(:, :, 2) + f(:, :, 6) + f(:, :, 9) ...
             - (f(:, :, 3) + f(:, :, 7) + f(:, :, 8));
    ux = rho_x_ux ./ rho;
    rho_x_uy = f(:, :, 4) + f(:, :, 6) + f(:, :, 7) ...
             - (f(:, :, 5) + f(:, :, 8) + f(:, :, 9));
    uy = rho_x_uy ./ rho;

    % Equilibrium distribution function for f
    for k = 1:ndir
        feq(:, :, k) = w(k)*rho;
    end

    
    % LBM for fluid temperature g second
    % Collision
    gcol = g - g * 1/tau_c + geq * 1/tau_c;
    
    % Streaming for g - Explicit version
    % This streaming implementation automatically applies periodic boundary
    % conditions in all edges of the computational domain.
    for k = 1:ndir
        for i = 1:ny
            for j = 1:nx
                xstreamed = mod(j + zeta_x(k), nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(i + zeta_y(k), ny);
                if ystreamed == 0
                    ystreamed = ny;
                end
                g(ystreamed, xstreamed, k) = gcol(i, j, k);
            end
        end
    end


    % Boundary conditions for g
    % LHS Boundary T = T_h using Inamuro
    for i = 1:ny
        for j = 1
            for k = [2 6 9]
                Tdash = (12/(2)) .* (T_h - g(i,j,1) - g(i,j,3) - g(i,j,4) ...
                                            - g(i,j,5) - g(i,j,7) - g(i,j,8));
                g(i,j,k) = w(k) * Tdash .* (1);
            end
        end
    end
    % RHS Boundary T = T_c using Inamuro
    for i = 1:ny
        for j = nx
            for k = [4 7 8]
                Tdash = (12/(2)) .* (T_c*1.0141 - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                         - g(i,j,5) - g(i,j,6) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash .* (1);
            end
        end
    end
    % % LHS Boundary T = T_h using Anti-BB
    % opp = [0 2 1 4 3 7 8 5 6];
    % for i = 1:ny
    %     for j = 1
    %         for k = [2 6 9]
    %             g(i,j,opp(k)+1) = -gcol(i,j,k) + 2*w(k) * T_h;
    %         end
    %     end
    % end
    % % RHS Boundary T = T_c using Anti-BB
    % for i = 1:ny
    %     for j = nx
    %         for k = [4 7 8]
    %             g(i,j,opp(k)+1) = -gcol(i,j,k) + 2*w(k) * T_c;
    %         end
    %     end
    % end
    % Top Boundary ie Anti-BB
    % opp = [0 2 1 4 3 7 8 5 6];
    % for i = ny
    %     for j = 1:nx
    %         for k = 1:ndir
    %             g(i,j,opp(k)+1) = -gcol(i,j,k) + 2*w(k) * C_p;
    %         end
    %     end
    % end
    % Bottom Boundary ie Anti-BB
    % opp = [0 2 1 4 3 7 8 5 6];
    % for i = 1
    %     for j = 1:nx
    %         for k = 1:ndir
    %             f(i,j,opp(k)+1) = -fcol(i,j,k) + 2*w(k) * C_p;
    %         end
    %     end
    % end


    % Macroscopic variables from g
    
    rho_x_T = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
            + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    % figure;
    % heatmap(T)
    T = rho_x_T./rho;
    % figure;
    % heatmap(T)


    % Equilibrium distribution function for g
    geq(:, :, 1) = 0;
    for k = 2:5
        geq(:, :, k) = rho.* T / 9 * 1.5;
    end
    for k = 6:9
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        geq(:, :, k) = rho.* T / 36 * 3;
    end


    % Force computation
    G(:,:,3) = rho(:,:) * beta * (T(:,:) - T_0).* grav;
    F(:,:,3) = (G(:,:,3).*zeta_x(3)./T) .* feq(:, :, 3);

    
    if mod(t, 100) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end


end


%% Saving T to file:
save /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_heat/T2 T









