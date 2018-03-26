clear;clc;
% 1 step
% 2 khist
% 3 nx
% 4 y
% 5 z
% 6 vx
% 7 vy
% 8 vz
% 9 Ep
% 10 Ek
% 11 T
% 12 Etot
% 13 type

% // MSE6270_MD uses the following units:
% // x - [A], t - [ps], m - [Da], E - [eV], F - [eV/A], P - [Pa]
% // define several constants for unit conversion
% static const double e = 1.60219e-19;       /*!< Electron charge */
% static const double kb = 8.617385e-5;      /*!< Boltzman constant [eV/K] */

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
    "equilT85", "T = 85 K"; ...
    "equilT90", "T = 90 K"; ...
    "equilT91", "T = 91 K"; ...
    "equilT92", "T = 92 K"; ...
    "equilT93", "T = 93 K"; ...
    "equilT94", "T = 94 K"; ...
    "equilT95", "T = 95 K"; ...
    "equilT100", "T = 100 K"; ...
%     "nptT60", "T = 60 K"; ...
%     "nptT74", "T = 74 K"; ...
%     "nptT75", "T = 75 K"; ...
%     "nptT76", "T = 76 K"; ...
%     "nptT77", "T = 77 K"; ...
%     "nptT78", "T = 78 K"; ...
%     "nptT79", "T = 79 K"; ...
%     "nptT80", "T = 80 K"; ...
%     "nptT85", "T = 85 K"; ...
%     "nptT90", "T = 90 K"; ...
% %     "nptT91", "T = 91 K"; ...
% %     "nptT92", "T = 92 K"; ...
% %     "nptT93", "T = 93 K"; ...
% %     "nptT94", "T = 94 K"; ...
%     "nptT95", "T = 95 K"; ...
%     "nptT100", "T = 100 K"; ...
];

datadir = "data";


% T = [];

for i = 1:length(directories)
    path = sprintf('%s/%s/', directories(i,1), datadir);
    syspath = sprintf('%s/%s.out', directories(i,1), directories(i,1));
    sys{i} = importdata(syspath,' ',17);
    files(i,:) = dir(sprintf('%s*.d', path));
    for k = 1:length(files(i,:))
        snapshot{i,k} = dlmread( ...
                        [path,files(i,k).name],' ');

        Ep(i,k) = sum(snapshot{i,k}(:,9));
        Ek(i,k) = sum(snapshot{i,k}(:,10));
        T(i,k) = sum(snapshot{i,k}(:,11));
        Etot(i,k) = sum(snapshot{i,k}(:,12));
        x1(i,k) = sum(snapshot{i,k}(:,11))/kb;
        speed{i,k} = aps2ms* sqrt( snapshot{i,k}(:,6).^2 + ...
                                   snapshot{i,k}(:,7).^2 + ...
                                   snapshot{i,k}(:,8).^2 );

    end
%     disp(k)
end

dt = 0.001; % psec
time = [0:5000*dt:100000*dt];

% %% Pressure vs Temperature
% figure; hold on; temps = []; pres = [];
% for i = 1:length(directories)
% % plot(sys{i}.data(:,2), ...
% %     [sys{i}.data(:,3) sys{i}.data(:,4) sys{i}.data(:,5)]);
%     temps = [temps; mean(sys{i}.data(:,6))];
%     pres = [pres; max(sys{i}.data(:,7)./(1e9))];
% end
% plot(temps, pres, 'LineWidth', 3); grid on;
% title('Temperature vs Pressure, final configuration');
% xlabel('Mean Temperature [K]','FontWeight','bold','Color','black');
% ylabel('Max Pressure [GPa]','FontSize',18,'FontWeight','bold','Color','black');
% xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% % saveas(gcf,'figures/q2-1.png');

% %% Temperature
% figure; hold on;
% itemp = []; ftemp = [];
% for i = 1:length(directories)
%     itemp = [itemp; sys{i}.data(1,6)];
%     ftemp = [ftemp; sys{i}.data(end,6)];
% % plot(sys{i}.data(:,2), ...
% %     [sys{i}.data(end,6)]);
% %     times = [times; mean(sys{i}.data(:,2))];
% %     temps = [temps; mean(sys{i}.data(:,6))];
%     
% %     pres = [pres; max(sys{i}.data(:,7)./(1e9))];
% end
% % plot(itemp,ftemp, 'LineWidth', 3); grid on;
% plot(itemp,ftemp, '*'); grid on;
% title('Temperature vs Pressure, final configuration');
% xlabel('Mean Temperature [K]','FontWeight','bold','Color','black');
% ylabel('Max Pressure [GPa]','FontSize',18,'FontWeight','bold','Color','black');
% xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% % saveas(gcf,'figures/q2-1.png');

% %% Total Energy Plot
% figure; plot(time, Ek+Ep, 'LineWidth', 3);
% legend(directories(:,2),'Location','NorthWest','Orientation','vertical');
% ylim([-30 -15]); grid on;
% title('Total Energy vs Time @ varying temperatures');
% xlabel('time [ps]','FontWeight','bold','Color','black');
% ylabel('Total Energy [eV]','FontSize',18,'FontWeight','bold','Color','black');
% xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% saveas(gcf,'figures/q2-1.png');

% %% Kinetic and potential energy plot
% figure; plot(time, Ek, time, Ep, 'LineWidth', 3);
% legend(directories(:,2),'Location','NorthWest','Orientation','vertical');
% % ylim([-30 -15]); grid on;
% title('Equilibrium Temperature @ varying initial temperatures');
% xlabel('time [ps]','FontWeight','bold','Color','black');
% ylabel('Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
% xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% % saveas(gcf,'figures/q2-1.png');

%% Equuilibration Temperature plot
figure; plot(time, 2*Ek/(3*N*kb), 'LineWidth', 3);
legend(directories(:,2),'Location','NorthWest','Orientation','vertical');
% ylim([-30 -15]); grid on;
title('Equilibrium Temperature @ varying initial temperatures');
xlabel('time [ps]','FontWeight','bold','Color','black');
ylabel('Temperature [K]','FontSize',18,'FontWeight','bold','Color','black');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% saveas(gcf,'figures/q2-1.png');

%%
% myT = 68;
% c = [0:.1:600];
% mbdist = 4*pi*(c.^2) .* ...
%           ((m_argon/(2*pi*kb_jk*myT))^(3/2)) .* ...
%           exp( (-m_argon*(c.^2)) / (2*kb_jk*myT) );
%       mbdist = 20*mbdist;
% 
%       
% % get velocities from last few timesteps      
% simvels = [speed{1,21}; speed{1,20}; speed{1,19}; speed{1,18}];
% 
% figure; hold on;
% hp = histogram(simvels,20,'Normalization','probability');
% plot(c,mbdist)  
      
      

























      
      
      
      
      
     
     