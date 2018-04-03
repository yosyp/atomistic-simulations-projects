%
% This corresponds to:
% HW5 Question 2
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
    for j = 1:N
%         if k > 1
%         if (snapshot{i,k}(j,3) - snapshot{i,k-1}(j,3)) > 20
%             pos_x(k,j) = snapshot{i,k}(j,3);
            
        pos_x(k,j) = snapshot{i,k}(j,3);
        pos_y(k,j) = snapshot{i,k}(j,4);
        pos_z(k,j) = snapshot{i,k}(j,5);
        if k == 1
            r0(j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
        else        
            r(k,j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
        end
    end
end

plot(pos_x(:,1))

%% Calculate Mean Square Displacement
msd(1) = 0;
for k = 2:length(files(i,:))
    for j = 1:N    
        msd(k) = msd(k-1) + ((r(k,j) - r0(j))^2)/N;
%         msd(k) = ((r(k,j) - r0(j))^2)/N;
    end
end
% disp(msd/6);
plot(msd)




