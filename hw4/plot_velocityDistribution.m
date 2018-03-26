clear;clc;

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
%     "nptT60", "T = 60 K"; ...
%     "nptT74", "T = 74 K"; ...
%     "nptT75", "T = 75 K"; ...
%     "nptT76", "T = 76 K"; ...
%     "nptT77", "T = 77 K"; ...
%     "nptT78", "T = 78 K"; ...
%     "nptT79", "T = 79 K"; ...
%     "nptT80", "T = 80 K", 80; ...
%     "nptT85", "T = 85 K", 85; ...
%     "nptT90", "T = 90 K", 90; ...
%     "nptT91", "T = 91 K", 91; ...
    "nptT92", "T = 92 K", 92; ...
%     "nptT93", "T = 93 K", 93; ...
%     "nptT94", "T = 94 K", 94;94 K ...
%     "nptT95", "T = 95 K", 95; ...
%     "nptT100", "T = 100 K", 100; ...
];

datadir = "data";
i = 1;
    path = sprintf('%s/%s/', directories(i,1), datadir);
    syspath = sprintf('%s/%s.out', directories(i,1), directories(i,1));
    sys{i} = importdata(syspath,' ',17);
    files(i,:) = dir(sprintf('%s*.d', path));
    for k = 1:length(files(i,:))
        snapshot{i,k} = dlmread( ...
                        [path,files(i,k).name],' ');
        speed(k,:) = aps2ms* sqrt( snapshot{i,k}(:,6).^2 + ...
                                   snapshot{i,k}(:,7).^2 + ...
                                   snapshot{i,k}(:,8).^2 );
        Ek(i,k) = sum(snapshot{i,k}(:,10));
    end

c = [0:.1:600];
N = 500;      


myT = mean(2*Ek/(3*N*kb));
mbdist = 4*pi*(c.^2) .* ...
          ((m_argon/(2*pi*kb_jk*myT))^(3/2)) .* ...
          exp( (-m_argon*(c.^2)) / (2*kb_jk*myT) );
mbdist = 12*mbdist;
simvels = speed(1:end,:);

figure; hold on;
hp = histogram(simvels,50,'Normalization','probability');
plot(c,mbdist, 'LineWidth', 8)  
    legend(sprintf('T = %4.2f K', myT),'Location','NorthWest','Orientation','vertical');
    title('Velocity Distribution vs Maxwell-Boltzmann Distribution');
    xlabel('Magnitude of Velocity [m/s]','FontWeight','bold','Color','black');
    ylabel('Probability [%/100]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
%     saveas(gcf,'figures/q4-90.png');




