%
% This corresponds to:
% HW5 Question 3
% 1 step    % 2 khist   % 3 x
% 4 y       % 5 z       % 6 vx
% 7 vy      % 8 vz      % 9 Ep
% 10 Ek     % 11 T      % 12 Etot
% 13 type
%
clear;clc;

% // MSE6270_MD uses the following units:
% // x - [A], t - [ps], m - [Da], E - [eV], F - [eV/A], P - [Pa]
% enunit = 1e4 * m0 / e;  % !< A constant used to convert [Da*A^2/ps] to [eV] 
kb = 8.617385e-5;       % !< Boltzman constant [eV/K] 
kb_jk = 1.38064852e-23; % !< Boltzman constant [J/K] 
m0 = 1.6605402e-27;     % !< Atomic mass unit [kg]
aps2ms = 100; % speed conversion: 1 A/ps = 100 m/s

m_argon = m0*39.948;    % units: kg
N = 1372;               % number of atoms 
box_size = 40.8045;     % size of periodic boundary condition box (in x,y,z)

directories = [
     "nptT138", "T = 138 K"; ...
%      "equilT138", "T = 138 K"; ...
];

i = 1;
datadir = "data";
path = sprintf('%s/%s/', directories(i,1), datadir);
files(1,:) = dir(sprintf('%s*.d', path));
for k = 1:length(files(i,:))
    snapshot{i,k} = dlmread( ...
                    [path,files(i,k).name],' ');
    for j = 1:N
        pos_x(k,j) = snapshot{i,k}(j,3);
        pos_y(k,j) = snapshot{i,k}(j,4);
        pos_z(k,j) = snapshot{i,k}(j,5);
        r(k,j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
    end
end

for j=1:N
    r_prev = r(1,j); % j-th atom first timestep = true coordinates   
    gt(1,j) = r(1,j);  % corrected cooridnates (1st timestep = true coordinates)
    % Loop through time trajectory of j-th atom
    for k = 2:length(pos_x(:,j))
        dr = r(k,j) - r_prev;       % displacement between current and previous timestep
        if dr > .5*box_size         % displacement too far "right"?
            dr = dr - box_size;     % replace displacement
        end
        if dr < -.5*box_size        % displacement too far "left"?
            dr = dr + box_size;     % replace displacement
        end
        gt(k,j) = gt(k-1,j) + dr;
        r_prev = r(k,j);
    end
end

% Plot Single Particle projection
for j=1:N
    prev_x = pos_x(1,j);
    prev_y = pos_y(1,j);
    prev_z = pos_z(1,j);
    new_x(1,j) = prev_x;
    new_y(1,j) = prev_y;
    new_z(1,j) = prev_z;
    for i=2:length(pos_x(:,j))
       dx = pos_x(i,j) - prev_x;
       if dx > .5*box_size
           dx = dx - box_size;
       elseif dx < -.5*box_size
           dx = dx + box_size;
       end
       new_x(i,j) = new_x(i-1,j) + dx;
       prev_x = pos_x(i,j);

       dy = pos_y(i,j) - prev_y;
       if dy > .5*box_size
           dy = dy - box_size;
       elseif dy < -.5*box_size
           dy = dy + box_size;
       end
       new_y(i,j) = new_y(i-1,j) + dy;
       prev_y = pos_y(i,j);   

       dz = pos_z(i,j) - prev_z;
       if dz > .5*box_size
           dz = dz - box_size;
       elseif dz < -.5*box_size
           dz = dz + box_size;
       end
       new_z(i,j) = new_z(i-1,j) + dz;
       prev_z = pos_z(i,j);     
    end
end

% Pick atom # to plot:
k = 123;

figure; hold on;
cmap = colormap; c = 1:length(new_x(:,k));
% change c into an index into the colormap
% min(c) -> 1, max(c) -> number of colors
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
% make a blank plot
plot3(new_x(:,k),new_y(:,k),new_z(:,k),'linestyle','none')
% add line segments
for j = 1:(length(new_x(:,k))-1)
    line(new_x(j:j+1,k),new_y(j:j+1,k),new_z(j:j+1,k),'color',cmap(c(j),:),'LineWidth', 3)
end
colorbar; grid on;
xlim([0 box_size]);
ylim([0 box_size]);
zlim([0 box_size]);
xlabel('x [A]','FontWeight','bold','Color','black');
ylabel('y [A]','FontWeight','bold','Color','black');
zlabel('z [A]','FontWeight','bold','Color','black');
grid on;  title('Trajectory Projection of Single Atom');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% saveas(gcf,'figures/q3-msd.png');

