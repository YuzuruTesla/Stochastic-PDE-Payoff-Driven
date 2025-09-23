clc; clear; close all;

L = 40;
maxt = 4e6;

x = linspace(0, L, 4001);
t = linspace(0, maxt, 2001);

% Parameters
V = 4; 
C = 6;

D_u_values = 4:0.05:9;

% Number of seeds
num_seeds = 50;

% Initialize arrays to store final p_H and p_D values for each D_u
average_p_H_over_seeds = zeros(size(D_u_values));
average_p_D_over_seeds = zeros(size(D_u_values));

% Loop over seeds
for seed = 1:num_seeds
    rng(seed, "twister");
    
    % Loop over D_u values
    for i = 1:length(D_u_values)
        D_u = D_u_values(i);
        m = 0; % Symmetry parameter

        sol = pdepe(m, @(x,t,u,DuDx)pdeFunction(x, t, u, DuDx, D_u), ...
                    @initialConditions, @boundaryConditions, x, t);

        % Extract solutions
        u_1 = sol(:,:,1);
        v_1 = sol(:,:,2);

        % Get the solution at t = maxt
        u_final = u_1(end, :);
        v_final = v_1(end, :);

        % Compute payoffs p_H and p_D at t = maxt
        frac_final = u_final ./ (u_final + v_final);
        frac_final(isnan(frac_final)) = 0; % Handle division by zero

        p_H_final = (V-C)/2 .* frac_final + V .* (1 - frac_final);
        p_D_final =      0 .* frac_final + V/2 .* (1 - frac_final);

        dx = x(2) - x(1);

        % Left Riemann sums (excluding the last element)
        total_p_H_left = sum(p_H_final(1:end-1) .* u_final(1:end-1)) * dx;
        total_p_D_left = sum(p_D_final(1:end-1) .* v_final(1:end-1)) * dx;
        total_u_left   = sum(u_final(1:end-1))             * dx;
        total_v_left   = sum(v_final(1:end-1))             * dx;

        % Right Riemann sums (excluding the first element)
        total_p_H_right = sum(p_H_final(2:end) .* u_final(2:end)) * dx;
        total_p_D_right = sum(p_D_final(2:end) .* v_final(2:end)) * dx;
        total_u_right   = sum(u_final(2:end))                 * dx;
        total_v_right   = sum(v_final(2:end))                 * dx;

        % Average Riemann sums
        total_p_H_avg = (total_p_H_left + total_p_H_right) / 2;
        total_p_D_avg = (total_p_D_left + total_p_D_right) / 2;
        total_u_avg   = (total_u_left   + total_u_right)   / 2;
        total_v_avg   = (total_v_left   + total_v_right)   / 2;

        % Compute the average payoffs
        average_p_H_final = total_p_H_avg / total_u_avg;
        average_p_D_final = total_p_D_avg / total_v_avg;

        % Accumulate average values over seeds
        average_p_H_over_seeds(i) = average_p_H_over_seeds(i) + average_p_H_final;
        average_p_D_over_seeds(i) = average_p_D_over_seeds(i) + average_p_D_final;
    end
end

% Compute the overall averages across all seeds
average_p_H_over_seeds = average_p_H_over_seeds / num_seeds;
average_p_D_over_seeds = average_p_D_over_seeds / num_seeds;

% Write to CSV: columns are [D_u, avg p_H, avg p_D]
csvwrite('NP_average_pHD_50.csv', [D_u_values' average_p_H_over_seeds' average_p_D_over_seeds']);


% PDE function 
function [c, f, s] = pdeFunction(x, t, u, DuDx, D_u)
    % Parameters
    D_v = 0.1;
    kappa = 0.001;
    V = 4;
    C = 6;
  
    % Solution variables
    u1 = u(1);
    u2 = u(2);
    
    % Spatial derivatives
    Du1Dx = DuDx(1);
    Du2Dx = DuDx(2);
    
    if u1 + u2 == 0
        frac = 0;
    else
        frac = u1 / (u1 + u2);
    end
    
    pi_u = (V-C)/2 * frac + V * (1 - frac);
    pi_v =      0 * frac + V/2 * (1 - frac);
    
    % PDE definitions
    c = [1; 1];
    f = [D_u * Du1Dx;
         D_v * Du2Dx];
    s = [u1 * (pi_u - kappa * (u1 + u2));
         u2 * (pi_v - kappa * (u1 + u2))];
end


function [pl, ql, pr, qr] = boundaryConditions(xl, ul, xr, ur, t)
    pl = [0; 0]; ql = [1; 1];  % zero‐flux left
    pr = [0; 0]; qr = [1; 1];  % zero‐flux right
end


function u0 = initialConditions(x)
    u0 = [4000/9 + 4*randn - 2; 2000/9 + 2*randn - 1];
end
