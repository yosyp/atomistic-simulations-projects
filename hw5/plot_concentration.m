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
% enunit = 1e4 * m0 / ex;  % !< A constant used to convert [Da*A^2/ps] to [eV] 
kb = 8.617385e-5;       % !< Boltzman constant [eV/K] 
kb_jk = 1.38064852e-23; % !< Boltzman constant [J/K] 
m0 = 1.6605402e-27;     % !< Atomic mass unit [kg]
aps2ms = 100; % speed conversion: 1 A/ps = 100 m/s

m_argon = m0*39.948;    % units: kg
N = 1280;               % number of atoms 

directories = [
     "part2", "T = 138 K"; ...
%      "equilpart2", "T = 138 K"; ...
];

i = 1;
% datadir = "data_gather2708";
datadir = "data";
path = sprintf('%s/%s/', directories(i,1), datadir);
files(1,:) = dir(sprintf('%s*.d', path));
for k = 1:length(files(i,:))
    snapshot{i,k} = dlmread( ...
                    [path,files(i,k).name],' ');
    for j = 1:N
        q(k,j) = snapshot{i,k}(j,2);
        x(k,j) = snapshot{i,k}(j,3);
        y(k,j) = snapshot{i,k}(j,4);
        z(k,j) = snapshot{i,k}(j,5);
    end
end

k=length(files(i,:));
% k=1;
% max = max( z(k,:) );
% % min = min( z(k,:) );
min = 0;
max = 115.6;
% bins = 11;
bins = 44;
dbin = (max-min)/(bins-1);
db = min:dbin:max;

% 11.56 
conc = zeros(length(db),3);
% conc = zeros(N,3);
for i=1:length(db)-1
    fprintf('more than %4.4f less than %4.4f\n', db(i), db(i)+dbin);
    for j=1:N
        if (z(k,j) >= db(i)) && (z(k,j) <= (db(i)+dbin))
            if q(k,j) == 1
                conc(i,1) = conc(i,1)+1;
            elseif q(k,j) == 3
                conc(i,2) = conc(i,2)+1;
            elseif q(k,j) == 0
                conc(i,3) = conc(i,3)+1;
            end
        end
    end
end

figure;
plot(db(1:end-1),conc(1:end-1,:),'LineWidth', 4);
title('Concentration Profile at t=961 ps');
legend('type A','RIGID','type B','location','North');
xlim([min max]);
xlabel('Z coordinate, [A]','FontWeight','bold','Color','black');
ylabel('# of atoms in plane','FontWeight','bold','Color','black');
grid on;
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);






% figure;plot(db,conc(:,1)+conc(:,3)); xlim([min max]);

