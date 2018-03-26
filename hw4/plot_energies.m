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
N = 500;                % number of atoms 

directories = [
%     "equilT60" ...
%     "equilT74" ...
%     "equilT75" ...
%     "equilT76" ...
%     "equilT77" ...
%     "equilT78" ...
%     "equilT79" ...
%     "equilT80" ...
    "equilT85", "T = 85 K", 85; ...
    "equilT90", "T = 90 K", 90; ...
    "equilT91", "T = 91 K", 91; ...
%     "equilT92", "T = 92 K", 92; ...
%     "equilT93", "T = 93 K", 93; ...
%     "equilT94", "T = 94 K", 94; ...
%     "equilT95", "T = 95 K", 95; ...
%     "equilT100", "T = 100 K", 100; ...
];

datadir = "data";
dt = 0.001; % psec
time = [0:5000*dt:100000*dt];
for i = 1:length(directories)
    path = sprintf('%s/%s/', directories(i,1), datadir);
    files(i,:) = dir(sprintf('%s*.d', path));
    for k = 1:length(files(i,:))
        snapshot{i,k} = dlmread( ...
                        [path,files(i,k).name],' ');
        T(i,k) = sum(snapshot{i,k}(:,11));                    
        Ek(i,k) = sum(snapshot{i,k}(:,10));
        Ep(i,k) = sum(snapshot{i,k}(:,9));
        Etot(i,k) = sum(snapshot{i,k}(:,12));

    end

end

% Equilibration Temperature plot
figure;
plot(time, [Ep; Ek; Etot; T/kb; Ek+Ep], 'LineWidth', 3);

    legend(directories(:,2),'Location','NorthEast','Orientation','vertical');
    title('Equilibrium Temperature @ varying initial temperatures');
    xlabel('time [ps]','FontWeight','bold','Color','black');
    ylabel('Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
%     saveas(gcf,'figures/q2-temp1.png');

% figure;
% plot(str2double(directories(:,3)), mean(2*Ek(:,end-4:end)'/(3*N*kb)), '-*', 'LineWidth', 3, 'MarkerSize', 12);
%     grid on; grid minor; title('Initial vs Equilibrated Temperatures');
%     xlabel('Initial Temperature [K]','FontWeight','bold','Color','black');
%     ylabel('Equilibrated Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
%     xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% %     saveas(gcf,'figures/q2-temp2.png');
% 
