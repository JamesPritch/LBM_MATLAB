%% Plotting of Model Finger with grid size 128
%% Read variables from file
load /Users/jpritch/Documents/MATLAB/model/130x500/saved/T3000
T = T(2:129,1:468); % This removes the upper and lower rows which weren't
                    % part of the actual model


%% Setting Variables
% Scalars
nx = 468;
ny = 128;
scale = 1.39e-1; % 18 is closest to mm in paper
a = 404; % Centre of curve at fingertip
T_c = -30; 
T_h = 26;
L_artery = round(1e-3/(1.78e-2/ny));

% Vectors
x = 0:nx-1;
y = (0:ny-1) - (ny-1)/2;


%% Cleaning up T
xs = zeros(ny,nx);
for i = 1:nx
    xs(:,i) = x(i)+1;
end
ys = zeros(ny,nx);
for i = 1:ny
    ys(i,:) = y(i);
end
r = zeros(ny,nx);
for i = 1:ny
    for j = 1:nx
        r(i,j) = sqrt((x(j)-a)^2+(y(i)+0.5)^2);
    end
end
[row,col] = find(abs(r)>(32*ny/64 - 0.5) & xs >= a);
for i = 1:length(row)
    T(row(i), col(i)) = NaN;
end

% Dimensionalising temperature
T = 75 * T + 230 - 273;

%% Plotting T
% Initial setup
% figure;
% axes('FontSize',11.5, 'NextPlot', 'add');
% % Plotting finger outline
% plot(zeros(ny),scale*((0:ny-1) - (ny-1)/2)/(ny-1), 'k', LineWidth=1);
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.5), 'k', LineWidth=1);
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.5), 'k', LineWidth=1);
% th = linspace( pi/2, -pi/2, 100);
% R = 1/2;
% x_c = R*cos(th) + (a)/(ny-1);
% y_c = R*sin(th) + 0;
% plot(scale*x_c,scale*y_c, 'k', LineWidth=1);
% % Plotting a
% plot(scale*(zeros(ny)+a/(ny-1)),scale*(2*(0:ny-1) - (ny-1))/(ny),...
%     'k:', LineWidth=1);
% % Plotting bone outline
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25), 'k', LineWidth=0.7);
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25), 'k', LineWidth=0.7);
% R = 1/4;
% x_c = R*cos(th) + (a)/(ny-1);
% y_c = R*sin(th) + 0;
% plot(scale*x_c,scale*y_c, 'k', LineWidth=0.7);
% % Plotting artery outline
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25 - L_artery/ny),...
%      'k', LineWidth=0.7);
% plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25 + L_artery/ny),...
%      'k', LineWidth=0.7);
% plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((ny/4-L_artery:ny/4)...
%     - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
% plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((3*ny/4-1:3*ny/4+L_artery-1)...
%     - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
% % Making plot look nice
% xlabel('${x/N_y}$','interpreter','latex', fontsize=19)
% ylabel('${y/N_y}$','interpreter','latex', fontsize=19)
% axis equal
% axis off
% xlim([-0.2*scale 3.85*scale])
% ylim([-0.7*scale 0.7*scale])

% Temperature heatmap
figure;
axes('FontSize',18, 'NextPlot', 'add');
% Plotting temperature of finger
s = pcolor(scale*x/(ny-1),scale*y/(ny-1),T);
set(s, 'EdgeColor', 'none');
% Plotting finger outline
plot(zeros(ny),scale*((0:ny-1) - (ny-1)/2)/(ny-1), 'k', LineWidth=1);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.5), 'k', LineWidth=1);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.5), 'k', LineWidth=1);
th = linspace( pi/2, -pi/2, 100);
R = 1/2;
x_c = R*cos(th) + (a)/(ny-1);
y_c = R*sin(th) + 0;
plot(scale*x_c,scale*y_c, 'k', LineWidth=1);
% Plotting bone outline
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25), 'k', LineWidth=0.7);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25), 'k', LineWidth=0.7);
R = 1/4;
x_c = R*cos(th) + (a)/(ny-1);
y_c = R*sin(th) + 0;
plot(scale*x_c,scale*y_c, 'k', LineWidth=0.7);
% Plotting artery outline
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25 - L_artery/ny),...
     'k', LineWidth=0.7);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25 + L_artery/ny),...
     'k', LineWidth=0.7);
plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((ny/4-L_artery:ny/4)...
    - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((3*ny/4-1:3*ny/4+L_artery-1)...
    - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
% Making plot look nice
xlabel('${x/N_y}$','interpreter','latex', fontsize=25)
ylabel('${y/N_y}$','interpreter','latex', fontsize=25)
axis equal
xlim([-0.2*scale 3.85*scale])
ylim([-0.7*scale 0.7*scale])
c = colorbar;
c.Limits = [T_c T_h];
box on

% Temperature contour
figure;
axes('FontSize',18, 'NextPlot', 'add');
% Plotting temperature of finger
contour(scale*x/(ny-1),scale*y/(ny-1),T, T_c:1:T_h, ...
                LineWidth=0.9);
% Plotting finger outline
plot(zeros(ny),scale*((0:ny-1) - (ny-1)/2)/(ny-1), 'k', LineWidth=1);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.5), 'k', LineWidth=1);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.5), 'k', LineWidth=1);
th = linspace( pi/2, -pi/2, 100);
R = 1/2;
x_c = R*cos(th) + (a)/(ny-1);
y_c = R*sin(th) + 0;
plot(scale*x_c,scale*y_c, 'k', LineWidth=1);
% Plotting bone outline
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25), 'k', LineWidth=0.7);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25), 'k', LineWidth=0.7);
R = 1/4;
x_c = R*cos(th) + (a)/(ny-1);
y_c = R*sin(th) + 0;
plot(scale*x_c,scale*y_c, 'k', LineWidth=0.7);
% Plotting artery outline
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) - 0.25 - L_artery/ny),...
     'k', LineWidth=0.7);
plot(scale*(0:a)/(ny-1),scale*(zeros(a+1) + 0.25 + L_artery/ny),...
     'k', LineWidth=0.7);
plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((ny/4-L_artery:ny/4)...
    - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
plot(scale*(zeros(L_artery+1)+a/(ny-1)),scale*((3*ny/4-1:3*ny/4+L_artery-1)...
    - (ny-1)/2)/(ny-1), 'k', LineWidth=0.7);
% Making plot look nice
xlabel('${x/N_y}$','interpreter','latex', fontsize=25)
ylabel('${y/N_y}$','interpreter','latex', fontsize=25)
axis equal
xlim([-0.2*scale 3.85*scale])
ylim([-0.7*scale 0.7*scale])
c = colorbar;
c.Limits = [T_c T_h];
box on











