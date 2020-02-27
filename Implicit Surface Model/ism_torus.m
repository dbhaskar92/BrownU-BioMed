%
% Agent-based model of particle movement on implicit surface (torus)
% Author: Dhananjay Bhaskar
% Last Modified: Jan 04, 2020
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
%

% number of particles
N = 150;

% params
Alpha = 1;
Sigma = 0.5;
phi = 1;
deltaT = 0.1;
totT = 20;

% init positions
X = zeros(N, 3);

% preallocate state variables
P = zeros(N, 3);
q = [2, 5];
Q = [-0.02, 0.025];
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 2);
dXdt = zeros(N, 3);
A = zeros(N, N);
b = zeros(N, 1);
lambda = zeros(N, 1);

% use rejection sampling for initial position
% https://math.stackexchange.com/questions/2017079/uniform-random-points-on-a-torus
cnt = 0;
r = q(1);
R = q(2);
while cnt < N
    U = rand();
    V = rand();
    W = rand();
    Theta = 2*pi*U;
    Phi = 2*pi*V;
    thresh = (R + r*cos(Theta))/(R + r);
    if W <= thresh
        cnt = cnt + 1;
        X(cnt, :) = [(R + r*cos(Theta))*cos(Phi), (R + r*cos(Theta))*sin(Phi), r*sin(Theta)]; 
    end
end

% visualize
visualize(X, r, R, 0);

t = 0;
itr = 0;
while t < totT

    % compute updated state vectors
    for i = 1 : N

        P(i,:) = [0, 0, 0]; 
        for j = 1 : N
            Fij = Alpha*exp(-1.0*norm((X(i,:)-X(j,:)))/(2*Sigma^2));
            P(i,:) = P(i,:) + (X(i,:) - X(j,:))*Fij;
        end

        F(i) = (X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)^2 - 4*q(2)^2*(X(i,1)^2 + X(i,2)^2);

        dFdX_i_x = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,1) - 8*q(2)^2*X(i,1);
        dFdX_i_y = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,2) - 8*q(2)^2*X(i,2);
        dFdX_i_z = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];

        dFdq_i_a = -4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(1);
        dFdq_i_R = 4*(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + q(2)^2 - q(1)^2)*q(2) - 8*q(2)*(X(i,1)^2 + X(i,2)^2); 
        dFdq(i,:) = [dFdq_i_a, dFdq_i_R];

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % assemble matrix A and vector b
    for i = 1 : N
        for j = 1 : N
            A(i,j) = dot(dFdq(i,:), dFdq(j,:));
        end
        b(i) = dot(dFdq(i,:), Q) + dot(dFdX(i,:), dXdt(i,:)) + phi*F(i);
    end
    
    % compute coefficients lambda
    lambda = A \ b;
    
    % compute rate of change of implicit surface
    weighted_sum = [0.0, 0.0];
    for j = 1 : N
        weighted_sum = weighted_sum + lambda(j)*dFdq(j,:);
    end
    dqdt = Q - weighted_sum;
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    % update geometry
    q(1) = q(1) + deltaT*dqdt(1);
    q(2) = q(2) + deltaT*dqdt(2);
    
    t = t + deltaT;
    itr = itr + 1;
    visualize(X, q(1), q(2), itr);
    
end

function [] = visualize(X, r, R, itr)

    theta_grid = linspace(0, 2*pi, 36);
    phi_grid = linspace(0, 2*pi, 18);
    [Phi_mesh, Theta_mesh] = meshgrid(phi_grid, theta_grid); 
    x = (R+r.*cos(Theta_mesh)).*cos(Phi_mesh);
    y = (R+r.*cos(Theta_mesh)).*sin(Phi_mesh);
    z = r.*sin(Theta_mesh);
    fig = figure('visible', 'off');
    mesh(x, y, z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim([-10 10])
    ylim([-8 8])
    zlim([-2 2])
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3));
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
    close

end