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

for j = 1:N % Loop through all atoms
    for i=2:length(pos_x(:,j))-1 % Loop through each timestep of 1 atom
% X        
        if abs(pos_x(i-1,j) - pos_x(i,j)) > box_size/2
            pos_x(i,j) = pos_x(i,j) - box_size;
        end        
        if abs(pos_x(i,j) - pos_x(i+1,j)) > box_size/2
            pos_x(i+1,j) = pos_x(i+1,j) - box_size;
        end       
        if abs(pos_x(i+1,j) - pos_x(i,j)) > box_size/2
            pos_x(i,j) = pos_x(i,j) + box_size;
        end

        
% Y        
        if abs(pos_y(i-1,j) - pos_y(i,j)) > box_size/2
            pos_y(i,j) = pos_y(i,j) + box_size;
        end    
        if abs(pos_y(i,j) - pos_y(i+1,j)) > box_size/2
            pos_y(i+1,j) = pos_y(i+1,j) + box_size;
        end     
        if abs(pos_y(i+1,j) - pos_y(i,j)) > box_size/2
            pos_y(i,j) = pos_y(i,j) + box_size;
        end     
       

% Z       
        if abs(pos_z(i-1,j) - pos_z(i,j)) > box_size/2
            pos_z(i,j) = pos_z(i,j) - box_size;
        end    
        if abs(pos_z(i,j) - pos_z(i+1,j)) > box_size/2
            pos_z(i+1,j) = pos_z(i+1,j) - box_size;
        end     
        if abs(pos_z(i+1,j) - pos_z(i,j)) > box_size/2
            pos_z(i,j) = pos_z(i,j) + box_size;
        end          
    end
end

% plot(pos_x(:,j)); ylim([-100 100]);

r0 = sqrt( (pos_x(1,:)).^2 + (pos_y(1,:)).^2 + (pos_z(1,:)).^2 );
r = sqrt( (pos_x(2:end,:)).^2 + (pos_y(2:end,:)).^2 + (pos_z(2:end,:)).^2 );



j = 32; figure;
subplot(2,2,1); plot(pos_x(:,j)); ylim([-100 100]); legend("x");
subplot(2,2,2); plot(pos_y(:,j)); ylim([-100 100]); legend("y");
subplot(2,2,3); plot(pos_z(:,j)); ylim([-100 100]); legend("z");
subplot(2,2,4); plot(r(:,j-1)); ylim([-100 100]); legend("r");


%% Calculate Mean Square Displacement
msd(1) = 0;
i=1;
for k = 2:length(files(i,:))-1
%     for j = 1:N    
%         msd(k) = msd(k-1) + ((r(k,j) - r0(j))^2)/N;
        msd(k) = sum( ((r(k,:) - r0).^2)/N );
%         fprintf("%d %d\n",k);
%     end
end
% disp(msd/6);

plot(msd/(6*100))

% Literature source for Argon
% at T = 84 K, 1.53 10e-5 cm^2/sec 
% https://aip.scitation.org/doi/abs/10.1063/1.1700899
%
% at T = 295 K and 42 kPa, 0.423 cm^2/sec
% https://journals.aps.org/pr/pdf/10.1103/PhysRev.72.1256

% Unit conversion:
% 1 Angstrom^2/ps = 10e-4 cm^2/sec = 10e-8 m^2/s






