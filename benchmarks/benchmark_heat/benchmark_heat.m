%% Lattice Boltzmann method code: convection in square cavity
%% Setting Grid Independant Variables
% Simulation parameters - input
% Scalars
nx = 64;
ny = 64;
niter = 50000;
T_0 = 0.5; % This is average between 0ºC and 20ºC ie 10ºC 
T_c = 0.45; % Set to non-dim and find proper non-dim later
T_h = 0.55;
rho_0 = 1.0000005;
beta = 1e-7; % Check hand calculation of this 100
grav = 1; % 9.8 in dim, idk for non-dim

% Vectors
x = 1:nx;
y = 1:ny;

% Matrices
rho = zeros(ny,nx) + rho_0; % Density

% D2Q9 velocity set parameters
ndir = 9;
c_s = sqrt(1/3);
c = sqrt(3*T_0); % Assuming T is constant throughout (it isn't)
zeta_x = c*[0, 1, -1, 0, 0, 1, -1, -1, 1];
zeta_y = c*[0, 0, 0, 1, -1, 1, 1, -1, -1];
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
udotu = ux.^2 + uy.^2;
for k = 1:ndir
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    feq(:, :, k) = w(k)*rho.*(1 + 3*zdotu/(c^2) + 9*zdotu.^2/(c^4) ...
                              - 1.5*udotu/(c^2));
end
f = feq;
fcol = zeros(ny, nx, ndir);

% Fluid temperature g
geq = zeros(ny, nx, ndir);
geq(:, :, 1) = 0;
for k = 2:5
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    geq(:, :, k) = rho.* T ./ 9 .* (1.5 + 1.5*zdotu./(c^2) + 2.25*zdotu.^2/(c^4) ...
                                  - 1.5*udotu./(c^2));
end
for k = 6:9
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    geq(:, :, k) = rho.* T ./ 36 .* (3 + 6*zdotu./(c^2) + 4.5*zdotu.^2/(c^4) ...
                                  - 1.5*udotu./(c^2));
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
                xstreamed = mod(j + zeta_x(k)/c, nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(i + zeta_y(k)/c, ny);
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
    %     f(i,1,:) = 0;
    % end
    % RHS Boundary ie u = 0
    % for i = 1:ny
    %     f(i,nx,:) = 0;
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
    udotu = ux.^2 + uy.^2;
    for k = 1:ndir
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        feq(:, :, k) = w(k)*rho.*(1 + 3*zdotu/(c^2) + 9*zdotu.^2/(c^4) ...
                                  - 1.5*udotu/(c^2));
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
                xstreamed = mod(j + zeta_x(k)/c, nx);
                if xstreamed == 0
                    xstreamed = nx;
                end
                ystreamed = mod(i + zeta_y(k)/c, ny);
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
                Tdash = (12/(2+3*uy(i))) .* (T_h - g(i,j,1) - g(i,j,3) - g(i,j,4) ...
                                                 - g(i,j,5) - g(i,j,7) - g(i,j,8));
                g(i,j,k) = w(k) * Tdash .* (1+3*uy(i));
            end
        end
    end
    % RHS Boundary T = T_c using Inamuro
    for i = 1:ny
        for j = nx
            for k = [4 7 8]
                Tdash = (12/(2+3*uy(i))) .* (T_c*1.0141 - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                            - g(i,j,5) - g(i,j,6) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash .* (1 + 3*uy(i));
            end
        end
    end
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
    udotu = ux.^2 + uy.^2;
    geq(:, :, 1) = -rho.* T * udotu / (3 * c^2);
    for k = 2:5
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        geq(:, :, k) = rho.* T / 9 .* (1.5 + 1.5*zdotu/(c^2) + 2.25*zdotu.^2/(c^4) ...
                                      - 1.5*udotu/(c^2));
    end
    for k = 6:9
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        geq(:, :, k) = rho.* T / 36 .* (3 + 6*zdotu/(c^2) + 4.5*zdotu.^2/(c^4) ...
                                      - 1.5*udotu/(c^2));
    end


    % Force computation
    G(:,:,3) = rho(:,:) * beta * (T(:,:) - T_0).* grav;
    F(:,:,3) = (G(:,:,3).*(zeta_x(3) - ux - uy)./T) .* feq(:, :, 3);

    
    if mod(t, 100) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end


end


%% Saving T to file:
save /Users/jpritch/Documents/MATLAB/benchmarks/benchmark_heat/T1 T









