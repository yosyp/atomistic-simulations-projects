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

directories = [
%     "equilT60" ...
%     "equilT74" ...
%     "equilT75" ...
%     "equilT76" ...
%     "equilT77" ...
%     "equilT78" ...
%     "equilT79" ...
%     "equilT80" ...
%     "equilT85" ...
    "nptT60" ...
%     "nptT74" ...
%     "nptT75" ...
%     "nptT76" ...
%     "nptT77" ...
%     "nptT78" ...
%     "nptT79" ...
%     "nptT80" ...
%     "nptT85" ...
];

datadir = "data";

T = [];

for i = 1:length(directories)
    path = sprintf('%s/%s/', directories(i), datadir);
    files(i,:) = dir(sprintf('%s*.d', path));
    for k = 1:length(files(i,:))
        snapshot{i,k} = dlmread( ...
                        [path,files(i,k).name],' ');

        Ep(i,k) = sum(snapshot{i,k}(:,9));
        Ek(i,k) = sum(snapshot{i,k}(:,10));
        T(i,k) = sum(snapshot{i,k}(:,11))/kb;
        speed{i,k} = aps2ms* sqrt( snapshot{i,k}(:,6).^2 + ...
                                   snapshot{i,k}(:,7).^2 + ...
                                   snapshot{i,k}(:,8).^2 );

    end
%     disp(k)
end


% 
% dt = 0.004;
% time = [0:dt:(length(Files)-1)*dt];
% 

myT = 68;
c = [0:.1:600];
mbdist = 4*pi*(c.^2) .* ...
          ((m_argon/(2*pi*kb_jk*myT))^(3/2)) .* ...
          exp( (-m_argon*(c.^2)) / (2*kb_jk*myT) );
      mbdist = 20*mbdist;

      
% get velocities from last few timesteps      
simvels = [speed{1,21}; speed{1,20}; speed{1,19}; speed{1,18}];

figure; hold on;
hp = histogram(simvels,20,'Normalization','probability');
plot(c,mbdist)  
      
      

























      
      
      
      
      
     
     