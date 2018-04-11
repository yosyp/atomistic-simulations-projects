%
% Make CNT's
% http://turin.nss.udel.edu/research/tubegenonline.html
% post-process:
% cut header
% $ awk '{ print $2 " " $3 " " $4}' input.txt  > output.txt

clear;clc;
% CNT = dlmread('cnt.xyz', ' ');
% CNT = dlmread('longcnt.xyz', ' ');
CNT = dlmread('midcnt.xyz', ' ');
x = CNT(:,1);
y = CNT(:,2);
z = CNT(:,3);
natoms = length(x);

% marker size
s = 250*ones(1,length(x))';

% subplot(2,2,1);
% hold on; grid on; view(40,35);
% xlim([-8 8]);
% ylim([-8 8]);
% zlim([-8 8]);
% scatter3(x,y,z,s,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]);

for i = 1:natoms
    for j = 1:natoms
        if i ~= j
            neighbs(j,i) = sqrt((x(j)-x(i))^2 + ... 
                              (y(j)-y(i))^2 + ... 
                              (z(j)-z(i))^2);
            neighbs(j,i) = round(neighbs(j,i), 4);
        else
            neighbs(j,i) = NaN;
        end
    end
end
% lowtriag = [tril(ones(natoms), 0), zeros(natoms)];
% neighbs(lowtriag == 1) = NaN;

bonds = zeros(natoms,natoms);
for i = 1:natoms
    nnindex = find(neighbs(i,:) == min(neighbs(i,:)));
    for k = nnindex
        line([x(i) x(k)], ...
             [y(i) y(k)], ...
             [z(i) z(k)], 'Color', 'red', 'LineWidth', 3);
         bonds(i,k) = 1;
         fprintf('Adding %d -> %d (%4.4f)\n', i , k, neighbs(i,k));         
    end
end
   
  
% subplot(2,2,2);
% image(bonds,'CDataMapping','scaled');
% 
% subplot(2,2,3);
% image(neighbs,'CDataMapping','scaled'); colorbar;

mygraph = graph(bonds);

firstnode = 24;
N1 = firstnode;
N2 = firstnode;
N3 = firstnode;

N1 = neighbors(mygraph,firstnode);
for i = 1:length(N1)
%     N2(i,:) = neighbors(mygraph, N1(i));
    N2 = [N2; neighbors(mygraph, N1(i))];
end
% N2 = reshape(N2,length(N2)^2,1);     % convert 3x3 to 9x1 vector
N2 = N2(N2 ~= firstnode);            % remove parent node from list

for i = 1:length(N2)
    nbs = neighbors(mygraph, N2(i));
    while length(nbs) < 3
        nbs = [nbs; NaN];
    end
    N3 = [N3; nbs];
end

% Do not backlink to first node
for i = 1:length(N1)
    N3 = N3(N3 ~= N1(i));
end
N3 = N3(N3 ~= firstnode);            % remove parent node from list

rings = zeros(natoms,natoms);
% for i = 1:length(N1)
i = 1;
    for j = 1:length(N2)
%         for k = 1:2:length(N3)
            if j == 1
                k = 1;
            else
                k = 2*j - 1;
            end

            my_ring = N1(i);            
            fprintf("i:%d j:%d j:%d, N1:%d N2:%d N3:%d\n", i,j,k,N1(i),N2(j),N3(k));
            rings(my_ring,firstnode) = 1;
            rings(my_ring,N1(i)) = 1;
            rings(my_ring,N2(j)) = 1;
            if isnan(N3(k+1))
                rings(my_ring,N3(k)) = 1;
            else
                dxk1 = sqrt((x(N3(k))-x(firstnode))^2 + ... 
                            (y(N3(k))-y(firstnode))^2 + ... 
                            (z(N3(k))-z(firstnode))^2);            
                dxk2 = sqrt((x(N3(k+1))-x(firstnode))^2 + ... 
                            (y(N3(k+1))-y(firstnode))^2 + ... 
                            (z(N3(k+1))-z(firstnode))^2);   
                fprintf("%d -> %d = %4.4f\n%d -> %d = %4.4f\n", N3(k),firstnode,dxk1,N3(k+1),firstnode,dxk2);
                if dxk1 > dxk2 
                    rings(my_ring,N3(k+1)) = 1;
                    fprintf("%d wins!\n\n", N3(k+1));
                else
                    rings(my_ring,N3(k)) = 1;
                    fprintf("%d wins!\n\n", N3(k));                    
                end
            end
            
%             disp(rings(my_ring,:));
%         end
    end
% end

image(rings,'CDataMapping','scaled');

% figure;plot(mygraph);
% [n,m] = size(N2);
% for i = 1:n
% %     k = 1;
%     for j = 1:m
% %         if N2(n,m) ~= firstnode
%         dx(i,j) = sqrt((x(firstnode)-x(N2(i,j)))^2 + ... 
%                        (y(firstnode)-y(N2(i,j)))^2 + ... 
%                        (z(firstnode)-z(N2(i,j)))^2);
% %             k = k+1;
% %         end
%     end
% end


% i = 20;
% firstn = find(neighbs(i,:) == min(neighbs(i,:)));
% 
% secondn = zeros(length(firstn),1);
% k = 1;
% for j = firstn
%     second(k,3) = find(neighbs(k,:) == min(neighbs(k,:)));
%     k = k + 1;
% end

% thirdn 

  