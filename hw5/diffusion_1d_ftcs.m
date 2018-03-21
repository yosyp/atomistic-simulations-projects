clear;clc;

D = .1;  % diffusion constant
L = 1;  % system size
dx = .05; % step of spatial discretization
h = .01;  % timestep
max_t = 10;  % total simulation time

N = L / dx; % Number of nodes
beta = (h*D/(dx^2)); % Constant used in main loop

if (2*D*h/(dx^2)) > .99
    disp('Von Neumann stability condition not met!');
end

% initial concentration profile
x = [0:dx:L];
sigma = 0.1; mu = .5;
C(1,:) = exp(-((x-mu).^2)/(2*sigma^2));

% loop over time
step = 2;
for t = 0:h:max_t
    for i = 2:N-1
        C(step,i) = C(step-1,i) + ...
                    beta*(C(step-1,i+1) - 2*C(step-1,i) + C(step-1,i-1));
    end
    C(step,1) = C(step-1,1) + 2*beta*(C(step-1,2) - C(step-1,1));
    C(step,N) = C(step-1,N) - 2*beta*(C(step-1,N) - C(step-1,N-1));
    step = step+1;
end

surf(C);