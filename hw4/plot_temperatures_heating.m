%
% This corresponds to:
% HW4 Question 2
%
clear;clc;
% 1 step    % 2 khist   % 3 nx
% 4 y       % 5 z       % 6 vx
% 7 vy      % 8 vz      % 9 Ep
% 10 Ek     % 11 T      % 12 Etot
% 13 type

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
%     "nptT60", "T = 60 K"; ...
%     "nptT74", "T = 74 K"; ...
%     "nptT75", "T = 75 K"; ...
%     "nptT76", "T = 76 K"; ...
%     "nptT77", "T = 77 K"; ...
%     "nptT78", "T = 78 K"; ...
%     "nptT79", "T = 79 K"; ...
    "nptT80", "T = 80 K", 80; ...
    "nptT85", "T = 85 K", 85; ...
    "nptT90", "T = 90 K", 90; ...
    "nptT91", "T = 91 K", 91; ...
    "nptT92", "T = 92 K", 92; ...
    "nptT93", "T = 93 K", 93; ...
    "nptT94", "T = 94 K", 94; ...
    "nptT95", "T = 95 K", 95; ...
%     "nptT100", "T = 100 K", 100; ...
];

datadir = "data";
dt = 0.001; % psec
time = [0:2500*dt:100000*dt];
for i = 1:length(directories)
    path = sprintf('%s/%s/', directories(i,1), datadir);
    files(i,:) = dir(sprintf('%s*.d', path));
    for k = 1:length(files(i,:))
        snapshot{i,k} = dlmread( ...
                        [path,files(i,k).name],' ');
        T(i,k) = sum(snapshot{i,k}(:,11));                    
        Ek(i,k) = sum(snapshot{i,k}(:,10));

    end

end

% Temperature plot
figure;
plot(time, 2*Ek/(3*N*kb), 'LineWidth', 3);

    legend(directories(:,2),'Location','NorthEast','Orientation','vertical');
    ylim([75 115]);
    title('Post-heating @ varying initial temperatures');
    xlabel('time [ps]','FontWeight','bold','Color','black');
    ylabel('Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
%     saveas(gcf,'figures/q2-tempheating1.png');

figure;
plot(str2double(directories(:,3)), mean(2*Ek(:,end-4:end)'/(3*N*kb)), '-*', 'LineWidth', 3, 'MarkerSize', 12);
    grid on; grid minor; title('Initial vs Post-Heating Temperatures');
    xlabel('Initial Temperature [K]','FontWeight','bold','Color','black');
    ylabel('Final Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
%     saveas(gcf,'figures/q2-tempheating2.png');

