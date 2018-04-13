clear;clc;


% T  = 1000;              % K
% T = [0.001 0.01 0.1 0 1 10 100 1000 1000];
% T = logspace(0, 4);
T = [10:1:1000];
kT = T*8.617385E-05;  % eV
N = 100*100;

evf = 0.05; % eV

neq = N*exp(-evf ./ kT);
semilogy(T, neq,'LineWidth', 6);

    title('Equilibrium Vacancy Concentration, \epsilon_v = 0.05 eV'); grid on;
    xlabel('Temperature [K]','FontWeight','bold','Color','black');
    ylabel('Equilibrium Concentration','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);


