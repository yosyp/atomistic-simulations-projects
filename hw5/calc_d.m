%
% This corresponds to:
% HW5 Question 2
% 1 step    % 2 khist   % 3 x
% 4 y       % 5 z       % 6 vx
% 7 vy      % 8 vz      % 9 Ep
% 10 Ek     % 11 T      % 12 Etot
% 13 type
%
% Reference slides:
% https://ocw.mit.edu/courses/materials-science-and-engineering/3-021j-introduction-to-modeling-and-simulation-spring-2012/part-i-lectures-readings/MIT3_021JS12_P1_L3.pdf
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
box_size = 40.8045;     % size of periodic boundary condition box (in x,y,z)

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
        pos_x(k,j) = snapshot{i,k}(j,3);
        pos_y(k,j) = snapshot{i,k}(j,4);
        pos_z(k,j) = snapshot{i,k}(j,5);
        r(k,j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
    end
end

% "Unwrap" coordinates from periodic boundary simulation
% This undoes the gather() command that keeps atoms in simulation box
% Displacement vector is unwrapped instead of individual coordinates
for j=1:N
    r_prev = r(1,j); % j-th atom first timestep = true coordinates   
    gt(1,j) = r(1,j);  % corrected cooridnates (1st timestep = true coordinates)
    % Loop through time trajectory of j-th atom
    for k = 2:length(pos_x(:,j))
        dr = r(k,j) - r_prev;       % displacement between current and previous timestep
        if dr > .5*box_size         % displacement too far "right"?
            dr = dr - box_size;     % replace displacement
        end
        if dr < -.5*box_size        % displacement too far "left"?
            dr = dr + box_size;     % replace displacement
        end
        gt(k,j) = gt(k-1,j) + dr;
        r_prev = r(k,j);
    end
end

% Calculate Mean Square Displacement
msd(1) = 0;
i=1; % All files from only this simulation run (see directories variable)
for k = 1:length(files(i,:))
        msd(k) = sum( ((gt(k,:) - gt(1,:)).^2)/N );
end
msd = msd./(6*100);

p = polyfit(1:100,msd(1:end),1);
f1 = polyval(p,1:100);

txt1 = '<\Delta r(t)^2> = A + 6Dt';
txt2 = sprintf('A = %4.4f\nD = %4.4f', p(2), p(1));

figure; hold on;
plot(0:99, msd, 'LineWidth', 9); ylim([0 15]);
plot(0:99, f1, 'LineWidth', 5); ylim([0 15]);
t1 = text(5,14,txt1); t1.FontSize = 20; t1.FontWeight = 'bold';
t2 = text(5,12,txt2); t2.FontSize = 20; t2.FontWeight = 'bold';
legend('MSD','Line Fit', 'Location', 'SouthEast');
xlabel('Time [picosecond]','FontWeight','bold','Color','black');
ylabel('Mean Squared Displacement [Angstrom^2]','FontWeight','bold','Color','black');
grid on; grid minor; title('Mean Squared Displacement');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);
% saveas(gcf,'figures/q2-msd.png');
   