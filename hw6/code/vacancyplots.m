clear;clc;
% close all;

filename = 'vEnNB1.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d1=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

filename = 'vEnNB100.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d100=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

filename = 'vEnNB1000.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d1000=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

step{1}  = d1(:,1);
Etot{1}  = d1(:,2);
bonds{1} = d1(:,3:5);

step{2}  = d100(:,1);
Etot{2}  = d100(:,2);
bonds{2} = d100(:,3:5);

step{3}  = d1000(:,1);
Etot{3}  = d1000(:,2);
bonds{3} = d1000(:,3:5);

figE = figure;
% figBonds = figure;

figure(figE);
    loglog(step{1},Etot{1}./100, ...
           step{2},Etot{2}./100, ...
           step{3},Etot{3}./100, ...
           'LineWidth', 6);
    title('Energy vs MC Moves at varying temperatures'); grid on;
    xlabel('Move steps [Moves]','FontWeight','bold','Color','black');
    ylabel('Energy [KeV]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
    legend('T = 1 K','T = 100 K','T = 1000 K','location','SouthWest');

% figure(figBonds);
% % i=2;
% loglog(step{1},bonds{1}(:,3), ...
%        step{2},bonds{2}(:,3), ...
%        step{3},bonds{3}(:,3), ...
%     'LineWidth', 6);
%     legend('T = 1 K','T = 100 K','T = 1000 K','location','SouthWest');
% %     legend('AA','BB','AB','location','NorthWest');
%     title('# of AB Bonds vs MC Moves (log-log)'); grid on;
%     xlabel('Move steps [Moves]','FontWeight','bold','Color','black');
%     ylabel('# of bonds [Count]','FontSize',18,'FontWeight','bold','Color','black');
%     xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% 


































