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
%         if k == 1
%             r0(j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
%         else        
            r(k,j) = sqrt( (pos_x(k,j))^2 + (pos_y(k,j))^2 + (pos_z(k,j))^2 );
%         end
    end
end

for j=1:N
% j=12;             % only use 2nd atom instead of all 1372 to test algorithm
% if r(1,j) > .5*box_size
%     r(1,j) = r(1,j) - .5*box_size;
% end
r_prev = r(1,j); % j-th atom first timestep = true coordinates   
gt(1,j) = r(1,j);  % corrected cooridnates (1st timestep = true coordinates)
% Loop through time trajectory of j-th atom
for k = 2:length(pos_x(:,j))
    dr = r(k,j) - r_prev;       % displacement between current and previous timestep
%     fprintf('step %d: displacement %4.4f\n',k,dr);    
    if dr > .5*box_size         % displacement too far "right"?
        dr = dr - box_size;     % replace displacement
%         fprintf('\t too big: %d: %4.4f corrected\n', k,dr);
    end
    if dr < -.5*box_size        % displacement too far "left"?
        dr = dr + box_size;     % replace displacement
%         fprintf('\t too small %d: %4.4f corrected\n', k,dr);            
    end
    gt(k,j) = gt(k-1,j) + dr;
    r_prev = r(k,j);
%     gt(k,j) = r_prev + dr;        % store new atom displacement as previous 
% %                                 atom displacement + corrected displacement
%     r_prev = gt(k,j);             % current timestep corrected displacement will
                                % be used as previous timestep in loop
end
end

% figure; plot(1:100, r(:,j), 1:100, gt); legend("original","corrected");

% for k=1:length(pos_x(:,j))-1
% %     fprintf('%d %4.4f\n', k, gt(k+1) - gt(k));
%    if gt(k+1) - gt(k) > .5*box_size
%        fprintf('wtf\n');
%    end
% end
% j = 123;
% % for j = 1:N % Loop through all atoms
%     prev_r = r0(j);
%     fixed_disp(1,j) = r0(j);
%     for i=2:length(pos_x(:,j)) % Loop through each timestep of 1 atom
%         dr = r(i,j) - prev_r;
%         new_disp = dr;
% %         new_disp = prev_r + dr;
%         if new_disp > .5*box_size
%             fprintf('%d: %4.4f\n', i, new_disp);
%             new_disp = new_disp - box_size;
%         elseif new_disp < -.5*box_size
%             fprintf('%d: %4.4f\n', i, new_disp);            
%             new_disp = new_disp + box_size;
%         end
%         fixed_disp(i,j) = prev_r + new_disp;
%         prev_r = fixed_disp(i,j);
%     end
% % end
% figure;
% plot(1:100, fixed_disp(:,j), 1:100, r(:,j)); legend('corrected','original');

% Calculate the displacement of particles in comparison with the previous snapshot, 
    
% and if this displacement is larger/smaller than 0.5*box_size/-0.5*box_size, 
% correct the displacement by adding/subtracting box_size 
% (assume that particle cannot move box_size/2 within 1 ps).
% 
% 
% Create an array for corrected coordinates and assign corrected coordinates 
% at the first time step to be equal to coordinates from the first snapshot. 
% 
% Following corrected coordinates are obtained by adding corresponding 
% displacements to the corrected coordinates on the previous time step. 


% % X        
%         if abs(pos_x(i-1,j) - pos_x(i,j)) > box_size/2
%             pos_x(i,j) = pos_x(i,j) - box_size;
%         end        
%         if abs(pos_x(i,j) - pos_x(i+1,j)) > box_size/2
%             pos_x(i+1,j) = pos_x(i+1,j) - box_size;
%         end       
%         if abs(pos_x(i+1,j) - pos_x(i,j)) > box_size/2
%             pos_x(i,j) = pos_x(i,j) + box_size;
%         end
% % Y        
%         if abs(pos_y(i-1,j) - pos_y(i,j)) > box_size/2
%             pos_y(i,j) = pos_y(i,j) + box_size;
%         end    
%         if abs(pos_y(i,j) - pos_y(i+1,j)) > box_size/2
%             pos_y(i+1,j) = pos_y(i+1,j) + box_size;
%         end     
%         if abs(pos_y(i+1,j) - pos_y(i,j)) > box_size/2
%             pos_y(i,j) = pos_y(i,j) + box_size;
%         end     
% % Z       
%         if abs(pos_z(i-1,j) - pos_z(i,j)) > box_size/2
%             pos_z(i,j) = pos_z(i,j) - box_size;
%         end    
%         if abs(pos_z(i,j) - pos_z(i+1,j)) > box_size/2
%             pos_z(i+1,j) = pos_z(i+1,j) - box_size;
%         end     
%         if abs(pos_z(i+1,j) - pos_z(i,j)) > box_size/2
%             pos_z(i,j) = pos_z(i,j) + box_size;
%         end          


% plot(pos_x(:,j)); ylim([-100 100]);

% r0 = sqrt( (pos_x(1,:)).^2 + (pos_y(1,:)).^2 + (pos_z(1,:)).^2 );
% r = sqrt( (pos_x(2:end,:)).^2 + (pos_y(2:end,:)).^2 + (pos_z(2:end,:)).^2 );



% j = 32; figure;
% subplot(2,2,1); plot(pos_x(:,j)); ylim([-100 100]); legend("x");
% subplot(2,2,2); plot(pos_y(:,j)); ylim([-100 100]); legend("y");
% subplot(2,2,3); plot(pos_z(:,j)); ylim([-100 100]); legend("z");
% subplot(2,2,4); plot(r(:,j-1)); ylim([-100 100]); legend("r");


%% Calculate Mean Square Displacement
msd(1) = 0;
i=1;
for k = 1:length(files(i,:))
        msd(k) = sum( ((gt(k,:) - gt(1,:)).^2)/N );
end
msd = msd./6;

p = polyfit(1:100,msd(1:end),1);
f1 = polyval(p,1:100);

txt1 = '<\Delta r(t)^2> = A + 6Dt';
txt2 = sprintf('A = %4.4f\nD = %4.4f', p(2), p(1));

figure; hold on;
plot(0:99, msd, 'LineWidth', 9); ylim([0 250]);
plot(0:99, f1, 'LineWidth', 5); ylim([0 250]);
t1 = text(5,200,txt1); t1.FontSize = 20; t1.FontWeight = 'bold';
t2 = text(5,170,txt2); t2.FontSize = 20; t2.FontWeight = 'bold';
legend('MSD','Line Fit', 'Location', 'SouthEast');
xlabel('Time [picosecond]','FontWeight','bold','Color','black');
ylabel('Mean Squared Displacement [Angstrom^2]','FontWeight','bold','Color','black');
    grid on; grid minor; title('Mean Squared Displacement');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);

    
    
% Literature source for Argon
% at T = 84 K, 1.53 10e-5 cm^2/sec 
% https://aip.scitation.org/doi/abs/10.1063/1.1700899
%
% at T = 295 K and 42 kPa, 0.423 cm^2/sec
% https://journals.aps.org/pr/pdf/10.1103/PhysRev.72.1256

% at T = K and 87 kPa, 2.07 10e-5 cm^2/sec

% Unit conversion:
% 1 Angstrom^2/ps = 10e-4 cm^2/sec = 10e-8 m^2/s


%% Plot Single Particle projection
for j=1:N
    prev_x = pos_x(1,j);
    prev_y = pos_y(1,j);
    prev_z = pos_z(1,j);
    new_x(1,j) = prev_x;
    new_y(1,j) = prev_y;
    new_z(1,j) = prev_z;
    for i=2:length(pos_x(:,j))
       dx = pos_x(i,j) - prev_x;
       if dx > .5*box_size
           dx = dx - box_size;
       elseif dx < -.5*box_size
           dx = dx + box_size;
       end
       new_x(i,j) = new_x(i-1,j) + dx;
       prev_x = pos_x(i,j);

       dy = pos_y(i,j) - prev_y;
       if dy > .5*box_size
           dy = dy - box_size;
       elseif dy < -.5*box_size
           dy = dy + box_size;
       end
       new_y(i,j) = new_y(i-1,j) + dy;
       prev_y = pos_y(i,j);   

       dz = pos_z(i,j) - prev_z;
       if dz > .5*box_size
           dz = dz - box_size;
       elseif dz < -.5*box_size
           dz = dz + box_size;
       end
       new_z(i,j) = new_z(i-1,j) + dz;
       prev_z = pos_z(i,j);     
    end
end

figure; hold on;

k = 123;
cmap = colormap; c = 1:length(new_x(:,k));
% change c into an index into the colormap
% min(c) -> 1, max(c) -> number of colors
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
% make a blank plot
plot3(new_x(:,k),new_y(:,k),new_z(:,k),'linestyle','none')
% add line segments
for j = 1:(length(new_x(:,k))-1)
    line(new_x(j:j+1,k),new_y(j:j+1,k),new_z(j:j+1,k),'color',cmap(c(j),:),'LineWidth', 3)
end
colorbar; grid on;
xlim([0 box_size]);
ylim([0 box_size]);
zlim([0 box_size]);
xlabel('x [A]','FontWeight','bold','Color','black');
ylabel('y [A]','FontWeight','bold','Color','black');
zlabel('z [A]','FontWeight','bold','Color','black');
    grid on;  title('Mean Squared Displacement');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);

k = 2;
cmap = colormap; c = 1:length(new_x(:,k));
% change c into an index into the colormap
% min(c) -> 1, max(c) -> number of colors
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
% make a blank plot
plot3(new_x(:,k),new_y(:,k),new_z(:,k),'linestyle','none')
% add line segments
for j = 1:(length(new_x(:,k))-1)
    line(new_x(j:j+1,k),new_y(j:j+1,k),new_z(j:j+1,k),'color',cmap(c(j),:),'LineWidth', 3)
end
colorbar; grid on;
xlim([0 box_size]);
ylim([0 box_size]);
zlim([0 box_size]);

k = 12;
cmap = colormap; c = 1:length(new_x(:,k));
% change c into an index into the colormap
% min(c) -> 1, max(c) -> number of colors
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
% make a blank plot
plot3(new_x(:,k),new_y(:,k),new_z(:,k),'linestyle','none')
% add line segments
for j = 1:(length(new_x(:,k))-1)
    line(new_x(j:j+1,k),new_y(j:j+1,k),new_z(j:j+1,k),'color',cmap(c(j),:),'LineWidth', 3)
end
colorbar; grid on;
xlim([0 box_size]);
ylim([0 box_size]);
zlim([0 box_size]);

% % Loop through time trajectory of j-th atom
% for k = 2:length(pos_x(:,j))
%     dr = r(k,j) - r_prev;       % displacement between current and previous timestep
% %     fprintf('step %d: displacement %4.4f\n',k,dr);    
%     if dr > .5*box_size         % displacement too far "right"?
%         dr = dr - box_size;     % replace displacement
% %         fprintf('\t too big: %d: %4.4f corrected\n', k,dr);
%     end
%     if dr < -.5*box_size        % displacement too far "left"?
%         dr = dr + box_size;     % replace displacement
% %         fprintf('\t too small %d: %4.4f corrected\n', k,dr);            
%     end
%     gt(k,j) = gt(k-1,j) + dr;
%     r_prev = r(k,j);
% %     gt(k,j) = r_prev + dr;        % store new atom displacement as previous 
% % %                                 atom displacement + corrected displacement
% %     r_prev = gt(k,j);             % current timestep corrected displacement will
%                                 % be used as previous timestep in loop
% end




