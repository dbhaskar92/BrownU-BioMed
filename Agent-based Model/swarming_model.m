% 
% Author: Dhananjay Bhaskar
% Last Modified: January 16, 2018
% Description: Collective cell migration and clustering
% References:
% Theodore Kolokolnikov (Dalhousie University)
% José A. Carrillo (Imperial College London)
%

% Number of particles
n = 50;
model_type = "Potential";

% Initial particle position
z = randn(1, n) + 1i * randn(1, n);

% Initial polarization vector
toggle_polarity = "on";
p = 0.25 * (1/n) * (randn(1, n) + 1i * randn(1, n));
p_offset = randi([0,2000], 1, n); 

% Mitosis params
toggle_cell_cycle = "on";
cell_cycle_duration = 50000;
cell_cycle_offset = randi([0,40000], 1, n);

toggle_periodic_bdy = "on";

t = 0;
dt = 0.05;

% Drawing interval
tnext = 0;

% Iteration number
itr = 1;

while (itr > 0)

    t = t + dt;
    dz = z * 0;

    for j = 1 : n

        k = [1:j-1, j+1:n];
        r = abs(z(j) - z(k));
        F = 0;

        % Explicit force based model
        if model_type == "Force"
            
            F = F + tanh((1-r) * 7.7) + 0.7;

        % Potential based model
        elseif model_type == "Potential"

            % Attraction-Repulsion Kernel
            cA = 0.32;
            cR = 0.40; 
            lA = 0.10;
            lR = 0.08;
            U = -cA*exp(-r./lA) + cR*exp(-r./lR);
            U_grad = (cA/lA)*exp(-r./lA) - (cR/lR)*exp(-r./lR);
            F = F - U_grad;
            
        end
        
        % Repolarization (force magnitude weighted by number of cells)
        if mod(p_offset(j), 2500) == 0
            p(j) = 0.25 * (1/n) * (randn() + 1i * randn());
        end
        
        % Compute velocity of particle
        dz(j) = 1/n * sum( F.*(z(j)-z(k))./r );
        if toggle_polarity == "on"
            dz(j) = dz(j) + p(j);
        end
        
        % Cell divison
        if mod(cell_cycle_offset(j), cell_cycle_duration) == 0

            % Set new position and velocity
            new_pos = randn(1, n+1);
            new_pos(1, 1:n) = z;
            new_pos(1, n+1) = z(j);
            new_vel = new_pos * 0;
            new_vel(1, 1:n) = dz;
            new_vel(1, n+1) = -dz(j);

            % Update cell cycle & polarization vectors
            new_p = randn(1, n+1);
            new_p(1, 1:n) = p;
            new_p(1, n+1) = -p(j);
            new_p_offset = randi([0,2000], 1, n+1);
            new_p_offset(1, 1:n) = p_offset;
            new_p_offset(1, j) = 0;
            new_p_offset(1, n+1) = 0;
            new_cycle_offset = randi([0,40000], 1, n+1);
            new_cycle_offset(1, 1:n) = cell_cycle_offset;
            new_cycle_offset(1, j) = 0;
            new_cycle_offset(1, n+1) = 0;

            % Update
            n = n + 1;
            z = new_pos;
            dz = new_vel;
            p = new_p;
            p_offset = new_p_offset;
            cell_cycle_offset = new_cycle_offset;
            
        end

    end

    z = z + dz * dt;
    
    % Periodic boundary conditions
    if toggle_periodic_bdy == "on"
        for j = 1 : n

            % X-axis
            if real(z(j)) > 2
                z(j) = - 2 + (real(z(j)) - 2) + 1i * imag(z(j));
            elseif real(z(j)) < -2
                z(j) = 2 - abs(-2 - real(z(j))) + 1i * imag(z(j));
            end

            % Y-axis
            if imag(z(j)) > 2
                z(j) = -2i + 1i * (imag(z(j)) - 2) + real(z(j));
            elseif imag(z(j)) < -2
                z(j) = 2i - 1i * (abs(-2 - imag(z(j)))) + real(z(j));
            end
             
        end
    end
    
    itr = itr + 1;
    p_offset = p_offset + 1;

    if toggle_cell_cycle == "on"
        cell_cycle_offset = cell_cycle_offset + 1;
    end

    % Plot
    if t > tnext

        tnext = tnext + dt * 100;

        plot(z, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', [1 .6 .6]);

        axis([-2 2 -2 2])
        % axis equal;

        title(sprintf('t=%g, radius=%g, #cells=%g', t, max(abs(mean(z)-z)), n));
        drawnow;

        %t_string = sprintf('%08d', itr-1);
        %saveas(gcf, strcat(t_string, '.png'));
        %pause;

    end
   
end
