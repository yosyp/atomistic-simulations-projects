clear;clc;
close all;

filename = 'EnNB_1.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d1=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

filename = 'EnNB100.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d100=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

filename = 'EnNB1000.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d1000=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
fclose(fid);

filename = 'EnNB10000.out';
fid = fopen(filename, 'rt');  %the 't' is important!
d10000=cell2mat(textscan(fid,'%f %f %f %f %f','collectoutput',1));
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

step{4}  = d10000(:,1);
Etot{4}  = d10000(:,2);
bonds{4} = d10000(:,3:5);


% figE = figure;
figBonds = figure;

% figure(figE);
%     loglog(step{1},Etot{1}./100, ...
%            step{2},Etot{2}./100, ...
%            step{3},Etot{3}./100, ...
%            step{4},Etot{4}./100, ...
%            'LineWidth', 6);
%     title('Energy vs MC Moves at varying temperatures'); grid on;
%     xlabel('Move steps [Moves]','FontWeight','bold','Color','black');
%     ylabel('Energy [KeV]','FontSize',18,'FontWeight','bold','Color','black');
%     xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
%     legend('T = 1 K','T = 100 K','T = 1000 K','T = 10000 K','location','SouthWest');

figure(figBonds);
% i=2;
loglog(step{1},bonds{1}, ...
    'LineWidth', 6);
%     legend('T = 1 K','T = 100 K','T = 1000 K','T = 10000 K','location','NorthWest');
    legend('AA','BB','AB','location','NorthWest');
    title('Q3: # of Bonds vs MC Moves, T = 1 K'); grid on;
    xlabel('Move steps [Moves]','FontWeight','bold','Color','black');
    ylabel('# of bonds [Count]','FontSize',18,'FontWeight','bold','Color','black');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);



































