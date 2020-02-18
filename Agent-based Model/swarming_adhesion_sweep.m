% 
% Author: Dhananjay Bhaskar
% Last Modified: Jun 12, 2018
% Description: Collective cell migration and clustering
% References:
% Theodore Kolokolnikov (Dalhousie University)
% James Carrillo (Imperial College London)
%

polarity_sweep = linspace(0.25,0.55,21).*1e-2;
adhesion_sweep = linspace(3,6,21).*1e-2;

conditions_adh = string(linspace(1,21,21));
conditions_pol = string(linspace(1,21,21));

for adh_vals = 1 : length(adhesion_sweep)
    
    for pol_vals = 1 : length(polarity_sweep)

        cond_str = strcat(conditions_adh(adh_vals), "_adh_", conditions_pol(pol_vals), "_pol");
        mkdir(char(cond_str));
        
        % Seed random number generator
        % rng(1007)
        rng shuffle
        
        % Params
        n = 100;                            % Number of particles
        boxsize = 10;                       % Simulation domain
        
        polarity_strength = polarity_sweep(pol_vals);   
        adh_strength = adhesion_sweep(adh_vals);  
        
        nbd_eps = 1.2;                      % Radius of neighborhood
        cell_cycle_duration = 50000;        % Duration of cell cycle
        polarity_duration = 2500;           % Time until repolarization
        contact_inhibition_threshold = 4;   % Cell density threshold to stop proliferation

        % Lennard-Jones attraction-repulsion kernel params
        cA = 4 * 0.5;
        cR = 4 * 2.0;       
        lA = 1.0;
        lR = 0.5;

        cell_cycle_offset = 0.8 * cell_cycle_duration;
        cell_polarity_offset = 0.8 * polarity_duration;

        % Toggles
        toggle_periodic_bdy = "on";
        toggle_polarity = "on";
        toggle_cell_cycle = "on";
        
        % Initialze
        z = get_init_pos(boxsize, n);                                                           % Particle position
        p = polarity_strength * (get_unif_rand(boxsize, n) + 1i * get_unif_rand(boxsize, n));   % Cell polarity 
        p_timer = randi([0, cell_polarity_offset], 1, n);                                       % Repolarization timer
        mitosis_timer = randi([0, cell_cycle_offset], 1, n);                                    % Division timer

        % Drawing controls
        t = 0;
        dt = 0.02;
        tnext = 0;

        % Iteration number
        itr = 0;
        end_time = 2000 * 1e2;

        while (itr <= end_time)

            t = t + dt;
            dz = z * 0;
            avg_num_neighbours = 0;
            avg_speed = 0;

            [num_nbd, xbordercells, xborder_unit_vec, xborder_radii] = find_neighbours(z, nbd_eps, -boxsize, boxsize, -boxsize, boxsize);

            for j = 1 : n

                k = [1:j-1, j+1:n];
                r = abs(z(j) - z(k));

                % Repolarization
                if mod(p_timer(j), polarity_duration) == 0
                    p(j) = polarity_strength * (get_unif_rand(boxsize, 1) + 1i * get_unif_rand(boxsize, 1));
                end

                % Calculate number of neighbours
                num_nbd_cells = num_nbd(j);
                avg_num_neighbours = avg_num_neighbours + num_nbd(j);
                nearest_neighbours = r < nbd_eps;
                neighbour_index = find(nearest_neighbours);

                dz(j) = 0;
                if toggle_polarity == "on"
                    dz(j) = dz(j) + p(j);
                end

                % Directed adhesion force and Lennard-Jones type potential
                adh_vec = 0;
                LJ_vec = 0;
                U = 0;
                for nn = 1 : length(neighbour_index)
                    idx = neighbour_index(nn);
                    ridx = idx;
                    if idx >= j
                        idx = idx + 1;
                    end
                    unit_vec = (z(idx) - z(j))/r(ridx);
                    U = U - cA * exp(-r(ridx)/lA) + cR * exp(-r(ridx)/lR);
                    F = (cA/lA) * exp(-r(ridx)/lA) - (cR/lR) * exp(-r(ridx)/lR);
                    LJ_vec = LJ_vec + F * unit_vec;
                    if r(ridx) > 1.0
                        adh_vec = adh_vec + adh_strength * unit_vec;
                    end
                end
                for nn = 1 : length(xbordercells)
                    idx = xbordercells(nn);
                    if j == idx
                        unit_vec = xborder_unit_vec(nn);
                        dist = xborder_radii(nn);
                        U = U - cA * exp(-dist/lA) + cR * exp(-dist/lR);
                        F = (cA/lA) * exp(-dist/lA) - (cR/lR) * exp(-dist/lR);
                        LJ_vec = LJ_vec + F * unit_vec;
                        if xborder_radii(nn) > 1.0
                            adh_vec = adh_vec + adh_strength * unit_vec;    
                        end
                    end
                end
                dz(j) = dz(j) + adh_vec + (1.0/100) * LJ_vec;

                avg_speed = avg_speed + abs(dz(j));

            end

            avg_num_neighbours = avg_num_neighbours/n;
            avg_speed = avg_speed/n;

            % Update position
            z = z + dz * dt;

            % Cell divison
            j = 1;
            while (j <= n)

                % Calculate number of neighbours
                num_nbd_cells = num_nbd(j);

                division_event = false;

                if (mod(mitosis_timer(j), cell_cycle_duration) == 0 && mitosis_timer(j) > 0)

                    % Contact inhibition of cell division
                    if num_nbd_cells >= contact_inhibition_threshold
                        
                        mitosis_timer(j) = 0;

                    else

                        % Set new position and velocity
                        new_pos = get_unif_rand(boxsize, n+1) + 1i * get_unif_rand(boxsize, n+1);
                        new_pos(1, 1:n) = z;
                        new_pos(1, n+1) = z(j) + get_unif_rand(0.2, 1) + 1i * get_unif_rand(0.2, 1);
                        new_vel = new_pos * 0;
                        new_vel(1, 1:n) = dz;
                        new_vel(1, n+1) = -dz(j);

                        % Update cell cycle & polarization vectors
                        new_p = polarity_strength * (get_unif_rand(boxsize, n+1) + 1i * get_unif_rand(boxsize, n+1));
                        new_p(1, 1:n) = p;
                        new_p(1, n+1) = -p(j);

                        new_p_timer = randi([0, cell_polarity_offset], 1, n+1);
                        new_p_timer(1, 1:n) = p_timer;
                        new_p_timer(1, j) = 0;
                        new_p_timer(1, n+1) = 0;

                        new_mitosis_timer = randi([0, cell_cycle_offset], 1, n+1);
                        new_mitosis_timer(1, 1:n) = mitosis_timer;
                        new_mitosis_timer(1, j) = 0;
                        new_mitosis_timer(1, n+1) = 0;

                        % Update
                        z = new_pos;
                        dz = new_vel;
                        p = new_p;
                        p_timer = new_p_timer;
                        mitosis_timer = new_mitosis_timer;

                        n = n + 1;
                        division_event = true;

                    end

                end

                if (division_event == true)
                    j = 1;
                    [num_nbd, xbordercells, xborder_unit_vec, xborder_radii] = find_neighbours(z, nbd_eps, -boxsize, boxsize, -boxsize, boxsize);
                else
                    j = j + 1;
                end 

            end
    
            % Periodic boundary conditions
            if toggle_periodic_bdy == "on"

                for j = 1 : n

                    % X-axis
                    if real(z(j)) > boxsize
                        z(j) = -1*boxsize + (real(z(j)) - boxsize) + 1i * imag(z(j));
                    elseif real(z(j)) < (-1*boxsize)
                        z(j) = boxsize - abs((-1*boxsize) - real(z(j))) + 1i * imag(z(j));
                    end

                    % Y-axis
                    if imag(z(j)) > boxsize
                        z(j) = (-1*boxsize)*1i + 1i * (imag(z(j)) - boxsize) + real(z(j));
                    elseif imag(z(j)) < (-1*boxsize)
                        z(j) = (boxsize)*1i - 1i * (abs((-1*boxsize) - imag(z(j)))) + real(z(j));
                    end

                end

            end

            p_timer = p_timer + 1;

            if toggle_cell_cycle == "on"
                mitosis_timer = mitosis_timer + 1;
            end
    
            [num_nbd, xbordercells, xborder_unit_vec, xborder_radii] = find_neighbours(z, nbd_eps, -boxsize, boxsize, -boxsize, boxsize);

            % Plot
            if itr == tnext

                tnext = tnext + 100;
                
                neighbor_count = [];
                edge_x_pair = [];
                edge_y_pair = [];
                
                for j = 1 : n
                    
                    neighbor_count = [neighbor_count; num_nbd(j)];
                    
                    k = [1:j-1, j+1:n];
                    r = abs(z(j) - z(k));
                    nearest_neighbours = r < nbd_eps;
                    neighbour_index = find(nearest_neighbours);

                    for nn = 1 : length(neighbour_index)
                        idx = neighbour_index(nn);
                        if idx >= j
                            idx = idx + 1;
                        end
                        x_vals = [real(z(j)) real(z(idx))];
                        y_vals = [imag(z(j)) imag(z(idx))];
                        edge_x_pair = [edge_x_pair; x_vals];
                        edge_y_pair = [edge_y_pair; y_vals];
                    end
                    
                end
                
                for j = 1 : length(xbordercells)
        
                    x = real(z(xbordercells(j)));
                    y = imag(z(xbordercells(j)));
                    u = real(xborder_unit_vec(j));
                    v = imag(xborder_unit_vec(j));
                    if u == 0
                        parametric_t = find_min_pos([(-boxsize - y)/v (boxsize - y)/v]);
                    elseif v == 0
                        parametric_t = find_min_pos([(-boxsize - x)/u (boxsize - x)/u]);
                    else
                        parametric_t = find_min_pos([(-boxsize - y)/v (boxsize - y)/v (-boxsize - x)/u (boxsize - x)/u]);
                    end
                    intersection_x = x + parametric_t*u;
                    intersection_y = y + parametric_t*v;
                    edge_x_pair = [edge_x_pair; [x intersection_x]];
                    edge_y_pair = [edge_y_pair; [y intersection_y]];

                end
                
                itr_string = sprintf('%06d', itr);
                
                csvfname = strcat(cond_str, filesep, 'Pos_', itr_string, '.dat');
                csvwrite(csvfname, z)
                csvfname = strcat(cond_str, filesep, 'Velocity_', itr_string, '.dat');
                csvwrite(csvfname, dz)
                
                csvfname = strcat(cond_str, filesep, 'EdgeX_', itr_string, '.dat');
                csvwrite(csvfname, edge_x_pair)
                csvfname = strcat(cond_str, filesep, 'EdgeY_', itr_string, '.dat');
                csvwrite(csvfname, edge_y_pair)
                
                csvfname = strcat(cond_str, filesep, 'Neighbors_', itr_string, '.dat');
                csvwrite(csvfname, neighbor_count)
                
                disp_nbh = avg_num_neighbours;
                disp_spd = avg_speed;
                st = sprintf('Frame T = %06d, Avg. # Neighbours = %0.3f, Avg. Speed = %0.3f, N = %d', itr, disp_nbh, disp_spd, n);
                disp(st);

            end
    
            itr = itr + 1;
   
        end
        
        close all;
        
    end
end

function [res] = get_unif_rand(limit, n)
    res = -limit + 2 * limit * rand(1, n);
end

function [res] = get_init_pos(limit, n)

    cnt = 0;
    res = [];

    while cnt < n

        % Pad by 0.7 so particles across periodic bdy are not too close
        r1 = -(limit-0.7) + 2 * (limit-0.7) * rand(1, 1);
        r2 = -(limit-0.7) + 2 * (limit-0.7) * rand(1, 1);
        new_particle_pos = r1 + 1i * r2;

        if cnt == 0
            res = [new_particle_pos];
            cnt = 1;
        else
            num_particles_inited = size(res, 2);
            too_close = 0;
            % Reject positions where particles are too close
            for j = 1 : num_particles_inited
                if abs(res(j) - new_particle_pos) < 0.9
                    too_close = 1;
                end
            end
            if too_close == 0
                res = [res new_particle_pos];
                cnt = cnt + 1;
            end
        end

    end

end

function [num_nbd, xbordercells, xborder_unit_vec, xborder_radii] = find_neighbours(z, epsl, xmin, xmax, ymin, ymax)

    N = size(z, 2);
    num_nbd = zeros(1, N);
    xbordercells = [];
    xborder_unit_vec = [];
    xborder_radii = [];

    for p = 1 : N
        
        k = [1:p-1, p+1:N];
        r = abs(z(p) - z(k));
        nearest_nbd = r < epsl;
        num_nbd(p) = num_nbd(p) + sum(nearest_nbd);

        x_p = 0;
        y_p = 0;
        if real(z(p)) - xmin < epsl
            x_p = real(z(p)) - xmin + xmax;
        end
        if imag(z(p)) - ymin < epsl
            y_p = imag(z(p)) - ymin + ymax;
        end
        if x_p == 0 && y_p == 0
            continue;
        elseif x_p > 0 && y_p == 0
            y_p = imag(z(p));
        elseif x_p == 0 && y_p > 0
            x_p = real(z(p));
        end
        phantom_pos = x_p + 1i * y_p;
        r = abs(phantom_pos - z(k));
        nearest_nbd_phantom = r < epsl;
        nearest_nbd_phantom_idx = find(nearest_nbd_phantom);
        num_nbd(p) = num_nbd(p) + sum(nearest_nbd_phantom);
        
        for q = 1 : length(nearest_nbd_phantom_idx)
            idx = nearest_nbd_phantom_idx(q);
            ridx = idx;
             if idx >= p
                idx = idx + 1;
             end
            num_nbd(idx) = num_nbd(idx) + 1;
            xbordercells = [xbordercells p idx];
            u_vec = (z(idx) - phantom_pos)/r(ridx);
            xborder_radii = [xborder_radii r(ridx) r(ridx)];
            xborder_unit_vec = [xborder_unit_vec u_vec -u_vec];
        end
        
    end

end

function [res] = find_min_pos(t_list)

    tmp = [];
    for p = 1 : length(t_list)
        if t_list(p) > 0
            tmp = [tmp t_list(p)];
        end
    end
    res = min(tmp);

end
