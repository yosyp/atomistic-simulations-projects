clear;clc;

D = .9371;        % diffusion constant
L =  115.6;       % system size
dx = 1.25;        % step of spatial discretization
h = 0.5;          % timestep
max_t = 961;      % total simulation time
x = [0:dx:L];     % x vector

N = ceil(L / dx);       % Number of nodes
beta = (h*D/(dx^2));    % Constant used in main loop

if (2*D*h/(dx^2)) > .99
    disp('Von Neumann stability condition not met!');
end

% Initial concentration profile
typeAlen = floor(length(x)/2);              % half-length of z-coord
Ca(1,:) = zeros(1,length(x));               % type A atoms initial
Ca(1,1:typeAlen+1) = 38*ones(1,typeAlen+1); % type A atoms initial

Cb(1,:) = zeros(1,length(x));                       % type B atoms initial
Cb(1,typeAlen:length(x)-1) = 38*ones(1,typeAlen+1); % type B atoms initial

% Loop over time
step = 2;
for t = 0:h:max_t
    for i = 2:N-1
        Ca(step,i) = Ca(step-1,i) + ...
                    beta*(Ca(step-1,i+1) - 2*Ca(step-1,i) + Ca(step-1,i-1));
        Cb(step,i) = Cb(step-1,i) + ...
                    beta*(Cb(step-1,i+1) - 2*Cb(step-1,i) + Cb(step-1,i-1));        
    end
    Ca(step,1) = Ca(step-1,1) + 2*beta*(Ca(step-1,2) - Ca(step-1,1));
    Cb(step,1) = Cb(step-1,1) + 2*beta*(Cb(step-1,2) - Cb(step-1,1));
    Ca(step,N) = Ca(step-1,N) - 2*beta*(Ca(step-1,N) - Ca(step-1,N-1));
    Cb(step,N) = Cb(step-1,N) - 2*beta*(Cb(step-1,N) - Cb(step-1,N-1));
    step = step+1;
end

figure; 
plot(x,Ca(end,:),x,Cb(end,:),'LineWidth',4)
title('Concentration Profile at t=961 ps');
legend('type A','type B','location','North');
xlabel('Z coordinate, [A]','FontWeight','bold','Color','black');
ylabel('# of atoms in plane','FontWeight','bold','Color','black');
grid on;
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);  set(gca, 'LineWidth', 2);


