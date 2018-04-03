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
    speed(k,:) = aps2ms* sqrt( snapshot{i,k}(:,6).^2 + ...
                               snapshot{i,k}(:,7).^2 + ...
                               snapshot{i,k}(:,8).^2 );
    Ek(i,k) = sum(snapshot{i,k}(:,10));
    for j = 1:N
        pos_x(k,j) = snapshot{i,k}(j,3);
        pos_y(k,j) = snapshot{i,k}(j,4);
        pos_z(k,j) = snapshot{i,k}(j,5);
    end
end

% i = 300; % number of atom to plot
% cutoff = 20;
i = 900; % number of atom to plot
cutoff = 0;
x = pos_x(1:end-cutoff,i)';
y = pos_y(1:end-cutoff,i)';
z = zeros(size(x)); col = x;
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
xlim([0 50]); ylim([0 50]); grid on; colorbar; caxis([0 50]);
title('Trajectory of 1 Particle, color = time');
xlabel('X Position [Angstrom]','FontWeight','bold','Color','black');
ylabel('Y Position [Angstrom]','FontWeight','bold','Color','black');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);


% plot(pos_x(1,i), pos_y(1,i),'-ro','MarkerSize', 12);
% plot(pos_x(2:end-cutoff,i), pos_y(2:end-cutoff,i),'-ko','MarkerSize', 12);
% plot(pos_x(end-cutoff,i), pos_y(end-cutoff,i),'-bo','MarkerSize', 12);


