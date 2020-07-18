%
% Agent-based model of particle movement on a genus 2 surface (double torus)
% Author: Dhananjay Bhaskar
% Last Modified: Jul 17, 2020
%
% Reference: Using Particles to Sample and Control Implicit Surfaces
% Andrew P. Witkin and Paul S. Heckbert, Proc. SIGGRAPH '94
%
% Additional References:
% https://math.stackexchange.com/questions/152256/implicit-equation-for-double-torus-genus-2-orientable-surface
% http://stanwagon.com/wagon/mathimages/htmllinks/mathimages_1.html
%
% Use ffmpeg to generate movie:
% ffmpeg -r 10 -i sim_%03d.png -vcodec libx264 -y -an sim_movie.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
%

% seed RNG
rng(1337)

% number of particles
N = 100;

% params
Alpha = 0.5;
Sigma = 0.2;
phi = 1;
deltaT = 0.1;
totT = 10;

% init positions
X = zeros(N, 3);

% preallocate state variables
P = zeros(N, 3);
q = 0.2;
Q = 0.0;
F = zeros(N, 1);
dFdX = zeros(N, 3);
dFdq = zeros(N, 1);
dXdt = zeros(N, 3);

% use rejection sampling for initial position
cnt = 0;
while cnt < N
    U = -1+2*rand();
    V = -1+2*rand();
    W = rand();
    det = -U^4 + 2*U^6 - U^8 + 2*U^2*V^2 - 2*U^4*V^2 - V^4 + q^2;
    if det > 0
        if W <= 0.5
            cnt = cnt + 1;
            X(cnt, :) = [U, V, sqrt(det)]; 
        else
            cnt = cnt + 1;
            X(cnt, :) = [U, V, -1*sqrt(det)];
        end
    end
end

% visualize IC
visualize_genus_2_mesh(X, q, 0);

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

        F(i) = X(i,1)^8 - 2*X(i,1)^6 + X(i,1)^4 - 2*X(i,1)^2*X(i,2)^2 + 2*X(i,1)^4*X(i,2)^2 + X(i,2)^4 + X(i,3)^2 -q^2;

        dFdX_i_x = 8*X(i,1)^7 - 12*X(i,1)^5 + 4*X(i,1)^3 -4*X(i,1)*X(i,2)^2 + 8*X(i,1)^3*X(i,2)^2;
        dFdX_i_y = -4*X(i,1)^2*X(i,2) + 4*X(i,1)^4*X(i,2) + 4*X(i,2)^3;
        dFdX_i_z = 2*X(i,3);
        dFdX(i,:) = [dFdX_i_x, dFdX_i_y, dFdX_i_z];
        
        dFdq(i) = -2*q;

        correction = (dot(dFdX(i,:), P(i,:)) + dot(dFdq(i,:), Q) + phi*F(i))/(norm(dFdX(i,:))^2);
        dXdt(i,:) = P(i,:) - correction*dFdX(i,:);

    end
    
    % update position
    for i = 1 : N
        X(i,:) = X(i,:) + deltaT*dXdt(i,:);
    end
    
    t = t + deltaT;
    itr = itr + 1;
    
    visualize_genus_2_mesh(X, q, itr);
    
end

function [] = visualize_genus_2_mesh(X, q, itr)

    fig = figure('visible', 'off');
    [xs, ys] = meshgrid(-20:.01:20);
    zs = sqrt(-1.*xs.^4 + 2.*xs.^6 - xs.^8 + 2.*xs.^2.*ys.^2 - 2.*xs.^4.*ys.^2 - ys.^4 + q^2);
    zs(~(zs==real(zs))) = nan;
    surf(xs, ys, zs, 'edgecolor', 'none', 'facecolor', 'interp');
    hold on;
    surf(xs, ys, -zs, 'edgecolor', 'none', 'facecolor', 'interp');
    alpha 0.5;
    [xs, ys] = meshgrid(-1:.06:1);
    zs = sqrt(-1.*xs.^4 + 2.*xs.^6 - xs.^8 + 2.*xs.^2.*ys.^2 - 2.*xs.^4.*ys.^2 - ys.^4 + q^2);
    zs(~(zs==real(zs))) = nan;
    surf(xs, ys, zs, 'edgecolor', 'k', 'facecolor', 'none', 'edgealpha', 1);
    surf(xs, ys, -zs, 'edgecolor', 'k', 'facecolor', 'none', 'edgealpha', 1);
    % daspect([1 1 1])
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3), 20, 'filled');
    xlim([-1.2 1.2])
    ylim([-1.2 1.2])
    zlim([-0.3 0.3])
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
    close

end