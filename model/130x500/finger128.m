%% Lattice Boltzmann method code: convection in square cavity
%% Setting Grid Independant Variables
% Simulation parameters - input
% Scalars
nx = 500;
ny = 130;
ny_s = 128; % Simulation ny
niter = 3000;
rho_0 = 1; % This is standard
T_0 = 0.5; % T<1 seems to work
beta = 0.000016; % Check hand calculation of this 100
a = 404; % Centre of curve at fingertip
L_artery = round(1e-3/(1.78e-2/ny));

% Boundary condition parameters
T_c = 0.175; 
T_h = 0.917;
T_end = T_h;
u_finger = 0;

% Vectors
x = 1:nx;
y = (0:ny-1) - ny/2;
xs = zeros(ny,nx);
for i = 1:nx
    xs(:,i) = x(i);
end
T_blood = zeros(a) + T_h;
T_fing = zeros(a) + T_h;

% Matrices
rho = zeros(ny,nx) + rho_0; % Density

% D2Q9 velocity set parameters
ndir = 9;
c_s = sqrt(3);
zeta_x = [0, 1, -1, 0, 0, 1, -1, -1, 1];
zeta_y = [0, 0, 0, 1, -1, 1, 1, -1, -1];
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Simulation parameters - built
tau_v = 1; % I don't know this, this is a guess, doesn't seem to make difference
tau_c_bone = 1.75;
tau_c_tissue = 0.641;

% Initialisation of fields at time t = 0
T = zeros(ny, nx) + T_h; % Temperature
ux = zeros(ny, nx); % Velocity in x direction
uy = zeros(ny, nx); % Velocity in y direction
omega = zeros(ny,nx); % Damage function

% Force fields
G = zeros(ny,nx,ndir);
F = zeros(ny,nx,ndir);

% Curved outlines
% Finger outline
mask = zeros(ny,nx);
mask(   2,1:a+1) = 1;
mask(ny-1,1:a+1) = 1;
r = zeros(ny,nx);
for i = 1:ny
    for j = 1:nx
        r(i,j) = sqrt((x(j)-a)^2+(y(i)+0.5)^2);
    end
end
[row1,col1] = find(abs(r)>31*ny_s/64 & xs >= a & abs(r)<31.9*ny_s/64);
for i = 1:length(row1)
    mask(row1(i), col1(i)) = 1;
end
for j = a+25*ny_s/64:a+31*ny_s/64
    for i = 11*ny_s/64:54*ny_s/64
        if mask(i,j) == 1 && mask(i,j+1) == 1
            mask(i,j+1) = 0;
        end
    end
end
% Bone outline
tau_c = zeros(ny,nx) + tau_c_tissue;
[row2,col2] = find(xs >= a & abs(r)<17*ny_s/64);
for i = 1:length(row2)
    tau_c(row2(i), col2(i)) = tau_c_bone;
end
tau_c(ny_s/4 + 1:3*ny_s/4 + 1,1:a) = tau_c_bone;
% Initialising Temperature field outside finger
[row3,col3] = find(abs(r)>31*ny_s/64 & xs >= a);
for i = 1:length(row3)
    T(row3(i), col3(i)) = T_c;
end

% Directions for Inamuro on a curved b.c.
dir = cell(length(row1));
for i = 1:length(row1)
    if row1(i) <= 0.5 * ny
        if mask(row1(i)  ,col1(i)-1) == 1 && mask(row1(i)  ,col1(i)+1) == 1
            dir{i} = {5;8;9};
        end
        if mask(row1(i)  ,col1(i)-1) == 1 && mask(row1(i)+1,col1(i)+1) == 1
            dir{i} = {5;8};
        end
        if mask(row1(i)-1,col1(i)-1) == 1 && mask(row1(i)+1,col1(i)+1) == 1
            dir{i} = {3;5;8};
        end
        if mask(row1(i)-1,col1(i)-1) == 1 && mask(row1(i)  ,col1(i)+1) == 1
            dir{i} = {3;5;8;9};
        end
        if mask(row1(i)-1,col1(i)-1) == 1 && mask(row1(i)+1,col1(i)  ) == 1
            dir{i} = {3;8};
        end
        if mask(row1(i)-1,col1(i)  ) == 1 && mask(row1(i)+1,col1(i)+1) == 1
            dir{i} = {3;5;7;8};
        end
        if mask(row1(i)-1,col1(i)  ) == 1 && mask(row1(i)+1,col1(i)  ) == 1
            dir{i} = {3;7;8};
        end
    end
    if row1(i) > 0.5 * ny
        if mask(row1(i)  ,col1(i)-1) == 1 && mask(row1(i)  ,col1(i)+1) == 1
            dir{i} = {4;6;7};
        end
        if mask(row1(i)+1,col1(i)-1) == 1 && mask(row1(i)  ,col1(i)+1) == 1
            dir{i} = {3;4;6;7};
        end
        if mask(row1(i)  ,col1(i)-1) == 1 && mask(row1(i)-1,col1(i)+1) == 1
            dir{i} = {4;7};
        end
        if mask(row1(i)+1,col1(i)-1) == 1 && mask(row1(i)-1,col1(i)+1) == 1
            dir{i} = {3;4;7};
        end
        if mask(row1(i)+1,col1(i)-1) == 1 && mask(row1(i)-1,col1(i)  ) == 1
            dir{i} = {3;7};
        end
        if mask(row1(i)+1,col1(i)  ) == 1 && mask(row1(i)-1,col1(i)+1) == 1
            dir{i} = {3;4;7;8};
        end
        if mask(row1(i)+1,col1(i)  ) == 1 && mask(row1(i)-1,col1(i)  ) == 1
            dir{i} = {3;7;8};
        end
    end
end

% Temperature boundary condition
B = (log(T_h)-log(T_c))/a;
for i = 1:nx
    T_fing(i) = exp(log(T_h) - B * i);
end


%% Simulating using LBM
% Initialisation of the particle distribution function
% Fluid density f
feq = zeros(ny, nx, ndir);
udotu = ux.^2 + uy.^2;
for k = 1:ndir
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    feq(:, :, k) = w(k)*rho.*(1 + 3  .*zdotu   /(c_s^2) ...
                                + 9  .*zdotu.^2/(c_s^4) ...
                                - 1.5.*udotu   /(c_s^2));
end
f = feq;
fcol = zeros(ny, nx, ndir);

% Fluid temperature g
geq = zeros(ny, nx, ndir);
geq(:, :, 1) = -rho.* T .* udotu / (3 * c_s^2);
for k = 2:5
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    geq(:, :, k) = rho.* T ./ 9 .* (1.5 + 1.5 .*zdotu   /(c_s^2) ...
                                        + 2.25.*zdotu.^2/(c_s^4) ...
                                        - 1.5 .*udotu   /(c_s^2));
end
for k = 6:9
    zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
    geq(:, :, k) = rho.* T ./ 36 .* (3 + 6  .*zdotu   /(c_s^2) ...
                                       + 4.5.*zdotu.^2/(c_s^4) ...
                                       - 1.5.*udotu   /(c_s^2));
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


    % Boundary conditions of tissue for f
    % LHS Boundary of tissue ie u = 0 using Zou & He
    rho_w = 1/(1-u_finger) * (f(:,1,1) + f(:,1,4) + f(:,1,5) + ...
                         2*(f(:,1,3) + f(:,1,7) + f(:,1,8)));
    f(:,1,2) = f(:,1,3) + 2/3 * rho_w * u_finger;
    f(:,1,6) = f(:,1,8) - 0.5 * (f(:,1,4) - f(:,1,5)) ...
               + 1/6 * rho_w * u_finger + 1/2 * rho_w * u_finger;
    f(:,1,9) = f(:,1,7) + 0.5 * (f(:,1,4) - f(:,1,5)) ...
               + 1/6 * rho_w * u_finger - 1/2 * rho_w * u_finger;

    % Tip of finger u = 0 using bounceback
    opp = [0 2 1 4 3 7 8 5 6];
    for i = 1:length(row1)
        for k = 1:ndir
            f(row1(i),col1(i),opp(k)+1) = fcol(row1(i),col1(i),k);
        end
    end
    
    % Note implementing ny at 2 and ny-1 since 1 and ny are empty
    % Top Boundary of finger ie u = 0 using Zou & He
    rho_n = 1/(1-u_finger) * (f(2,:,1) + f(2,:,2) + f(2,:,3) + ...
                         2*(f(2,:,4) + f(2,:,6) + f(2,:,7)));
    f(2,:,5) = f(2,:,4) - 2/3 * rho_n * u_finger;
    f(2,:,8) = f(2,:,6) + 0.5 * (f(2,:,2) - f(2,:,3)) ...
               - 1/2 * rho_n * u_finger - 1/6 * rho_n * u_finger;
    f(2,:,9) = f(2,:,7) - 0.5 * (f(2,:,2) - f(2,:,3)) ...
               + 1/2 * rho_n * u_finger - 1/6 * rho_n * u_finger;
    
    % Bottom Boundary of finger ie u = 0 using Zou & He
    rho_s = 1/(1-u_finger) * (f(ny-1,:,1) + f(ny-1,:,2) + f(ny-1,:,3) + ...
                         2*(f(ny-1,:,5) + f(ny-1,:,8) + f(ny-1,:,9)));
    f(ny-1,:,4) = f(ny-1,:,5) + 2/3 * rho_s * u_finger;
    f(ny-1,:,6) = f(ny-1,:,8) - 0.5 * (f(ny-1,:,2) - f(ny-1,:,3)) ...
               + 1/2 * rho_s * u_finger + 1/6 * rho_s * u_finger;
    f(ny-1,:,7) = f(ny-1,:,9) + 0.5 * (f(ny-1,:,2) - f(ny-1,:,3)) ...
               - 1/2 * rho_s * u_finger + 1/6 * rho_s * u_finger;
    
    % Top blood vessel walls are u = 0 using bounceback
    opp = [0 2 1 4 3 7 8 5 6];
    for i = [3*ny_s/4 + 1 3*ny_s/4 + L_artery + 1]
        for j = 2:a
            for k = 1:ndir
                f(i,j,opp(k)+1) = fcol(i,j,k);
            end
        end
    end

    % Bottom blood vessel walls are u = 0 using bounceback
    for i = [ny_s/4 + 1 ny_s/4 + L_artery + 1]
        for j = 2:a
            for k = 1:ndir
                f(i,j,opp(k)+1) = fcol(i,j,k);
            end
        end
    end


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
        feq(:, :, k) = w(k)*rho.*(1 + 3  .*zdotu   /(c_s^2) ...
                                    + 9  .*zdotu.^2/(c_s^4) ...
                                    - 1.5.*udotu   /(c_s^2));
    end

    

    % LBM for fluid temperature g second
    % Collision
    gcol = g - g * 1./tau_c + geq * 1./tau_c;
    

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
    % LHS domain Boundary T = T_h using Inamuro
    for i = 1:ny
        for j = 1
            for k = [2 6 9]
                Tdash = (12/(2+3*uy(i))) .* (T_h - g(i,j,1) - g(i,j,3) - g(i,j,4) ...
                                                 - g(i,j,5) - g(i,j,7) - g(i,j,8));
                g(i,j,k) = w(k) * Tdash .* (1+3*uy(i));
            end
        end
    end
    
    % RHS domain Boundary T = T_c using Inamuro note u = 0 automatically
    for i = 1:ny
        for j = nx-5 % Shifted a couple to left for numerical stability
            for k = [4 7 8]
                Tdash = (12/(2+3*uy(i))) .* (T_c - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                 - g(i,j,5) - g(i,j,6) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % RHS finger curved Boundary T = T_c using Inamuro
    for i = 1:length(col1)
        T_in = 0;
        for j = 1:length(dir{i})
            T_in = T_in + g(row1(i),col1(i),dir{i}{j});
        end
        T_tot = g(row1(i),col1(i), 1) + g(row1(i),col1(i), 2) + g(row1(i),col1(i), 3) ...
              + g(row1(i),col1(i), 4) + g(row1(i),col1(i), 5) + g(row1(i),col1(i), 6) ...
              + g(row1(i),col1(i), 7) + g(row1(i),col1(i), 8) + g(row1(i),col1(i), 9);
        weight = 0;
        for j = 1:length(dir{i})
            weight = weight + w(dir{i}{j});
        end
        Tdash = (T_c - T_tot + T_in) / weight; % this temp doesn't work
        for j = cell2mat(dir{i})
            g(row1(i),col1(i),j) = w(j) * Tdash;
        end
    end

    % Top domain Boundary T = T_c (varying) using Inamuro. Note zero velocity
    for i = 2
        for j = 2:nx
            for k = [5 8 9]
                Tdash = 6 * (T_fing(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                 - g(i,j,4) - g(i,j,6) - g(i,j,7));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % Bottom domain Boundary T = T_c (varying) using Inamuro. Note zero velocity
    for i = ny-1
        for j = 2:nx
            for k = [4 6 7]
                Tdash = 6 * (T_fing(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                 - g(i,j,5) - g(i,j,8) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % Blood vessel boundary temperature
    lambda = (T_h - T_end)/a;
    for i = 1:a
        T_blood(i) = T_h - i * lambda;
    end

    % Upper blood vessel top wall (enforcing T upwards)
    for i = 3*ny_s/4 + L_artery + 1
        for j = 2:(a)
            for k = [4 6 7]
                Tdash = 6 * (T_blood(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                        - g(i,j,5) - g(i,j,8) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % Upper blood vessel bottom wall (enforcing T downwards)
    for i = 3*ny_s/4 + 1
        for j = 2:(a)
            for k = [5 8 9]
                Tdash = 6 * (T_blood(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                        - g(i,j,4) - g(i,j,6) - g(i,j,7));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % Lower blood vessel top wall (enforcing T upwards)
    for i = ny_s/4+1
        for j = 2:(a)
            for k = [4 6 7]
                Tdash = 6 * (T_blood(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                        - g(i,j,5) - g(i,j,8) - g(i,j,9));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end

    % Lower blood vessel bottom wall (enforcing T downwards)
    for i = ny_s/4 - L_artery + 1
        for j = 2:(a)
            for k = [5 8 9]
                Tdash = 6 * (T_blood(j) - g(i,j,1) - g(i,j,2) - g(i,j,3) ...
                                        - g(i,j,4) - g(i,j,6) - g(i,j,7));
                g(i,j,k) = w(k) * Tdash;
            end
        end
    end


    % Macroscopic variables from g
    rho_x_T = g(:, :, 1) + g(:, :, 2) + g(:, :, 3) + g(:, :, 4) + g(:, :, 5)...
            + g(:, :, 6) + g(:, :, 7) + g(:, :, 8) + g(:, :, 9);
    T = rho_x_T./rho;


    % Update sample temperature at boundary of bone to soft tissue
    T_end = T(ny_s/2,(a + 17*ny_s/64));
    T_end = T_end * (1+(T_h-T_end)/8.75);
    

    % Equilibrium distribution function for g
    udotu = ux.^2 + uy.^2;
    geq(:, :, 1) = -rho.* T .* udotu / (3 * c_s^2);
    for k = 2:5
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        geq(:, :, k) = rho.* T ./ 9 .* (1.5 + 1.5 .*zdotu   /(c_s^2) ...
                                            + 2.25.*zdotu.^2/(c_s^4) ...
                                            - 1.5 .*udotu   /(c_s^2));
    end
    for k = 6:9
        zdotu = zeta_x(k)*ux + zeta_y(k)*uy;
        geq(:, :, k) = rho.* T ./ 36 .* (3 + 6  .*zdotu   /(c_s^2) ...
                                           + 4.5.*zdotu.^2/(c_s^4) ...
                                           - 1.5.*udotu   /(c_s^2));
    end


    % Force computation
    G(:,:,3) = rho(:,:) * beta .* (T(:,:) - T_0);
    F(:,:,3) = (G(:,:,3).*(zeta_x(3) - ux - uy)./T) .* feq(:, :, 3);


    % Print progress
    if mod(t, 100) == 0
        fprintf('Iteration: %d, Time: %f \n', t, toc);
    end
    

    % Save variables to file
    if mod(t, 100) == 0
        T_name = sprintf('/Users/jpritch/Documents/MATLAB/model/130x500/saved/T%04d', t);
        save(T_name, 'T');
    end


end













