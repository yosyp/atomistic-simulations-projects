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

j = 2;
% my_x(1) = pos_x(1,1);
figure; hold on; plot(pos_x(:,j));
% my_x = pos_x(:,233);
% my_y = pos_y(:,1);
% my_z = pos_z(:,1);
box_size = 40.8045;
% for j = 1:N
    for i=2:length(pos_x(:,j))-1
        
        fprintf("%d: %4.4f\n", i, (pos_x(i+1,j) - pos_x(i,j)));
        if (pos_x(i+1,j) - pos_x(i,j)) > box_size/2
            fprintf("\telse %d: %4.4f\n", i, (pos_x(i+1,j) - pos_x(i,j)));
            pos_x(i,j) = pos_x(i,j) + box_size;
        end
        
        fprintf("%d: %4.4f\n", i, (pos_x(i,j) - pos_x(i+1,j)));        
        if (pos_x(i,j) - pos_x(i+1,:)) > box_size/2
            pos_x(i+1,j) = pos_x(i+1,j) + box_size;
            fprintf("\tif %d: %4.4f\n", i, (pos_x(i,j) - pos_x(i+1,j)));
        end
        
        fprintf("%d: %4.4f\n", i, (pos_x(i-1,j) - pos_x(i,j)));
        if (pos_x(i-1,j) - pos_x(i,j)) > box_size/2
            fprintf("\telse %d: %4.4f\n", i, (pos_x(i-1,j) - pos_x(i,j)));
            pos_x(i,j) = pos_x(i,j) + box_size;
        end        

%         if abs(pos_y(i,j) - pos_y(i+1,j)) > box_size/2
%             if ( abs(pos_y(i,j)) - abs(pos_y(i+1,j))) > (abs(pos_y(i+1,j)) - abs(pos_y(i,j)))
%                 pos_y(i+1,j) = pos_y(i+1,j) + box_size;                
%             else
%                 pos_y(i,j) = pos_y(i,j) + box_size; 
%             end
%         end     
%         if abs(pos_z(i,j) - pos_z(i+1,j)) > box_size/2
%             if ( abs(pos_z(i,j)) - abs(pos_z(i+1,j))) > (abs(pos_z(i+1,j)) - abs(pos_z(i,j)))
%                 pos_z(i+1,j) = pos_z(i+1,j) + box_size;                
%             else
%                 pos_z(i,j) = pos_z(i,j) + box_size; 
%             end
%         end        
    end
% end
plot(pos_x(:,j)); ylim([0 100]);

% r0 = sqrt( (pos_x(1,:)).^2 + (pos_y(1,:)).^2 + (pos_z(1,:)).^2 );
% r = sqrt( (pos_x(2:end,:)).^2 + (pos_y(2:end,:)).^2 + (pos_z(2:end,:)).^2 );
% figure; hold on;
% plot(pos_x(:,2));
% plot(pos_y(:,2));
% plot(pos_z(:,2));
% figure;plot(r(:,1));


% figure; hold on; plot(pos_x(:,1)); plot(my_x);
% figure; hold on; plot(pos_y(:,1)); plot(my_y);
% figure; hold on; plot(pos_z(:,1)); plot(my_z);

%% Calculate Mean Square Displacement
msd(1) = 0;
for k = 2:length(files(i,:))
    for j = 1:N    
%         msd(k) = msd(k-1) + ((r(k,j) - r0(j))^2)/N;
        msd(k) = ((r(k,j) - r0(j))^2)/N;
    end
end
% disp(msd/6);
plot(-msd)




