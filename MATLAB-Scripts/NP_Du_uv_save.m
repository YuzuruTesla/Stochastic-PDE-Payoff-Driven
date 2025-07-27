clc; clear; close all;

L = 40;
maxt = 4e6;

x = linspace(0, L, 4001);
t = linspace(0, maxt, 201);

% Parameters
V = 4; 
C = 6;

D_u_values = 4:0.05:9;

% Number of seeds
num_seeds = 50;

% Initialize arrays to store final u and v values for each seed and chi_u
average_u_final_over_seeds = zeros(size(D_u_values));
average_v_final_over_seeds = zeros(size(D_u_values));

% To store all u_final and v_final values for each seed and chi_u
u_final_all_seeds = zeros(num_seeds, length(D_u_values), length(x));
v_final_all_seeds = zeros(num_seeds, length(D_u_values), length(x));

% Loop over seeds
for seed = 1:num_seeds
    
    rng(seed, "twister");
    
    % Loop over chi_u values
    for i = 1:length(D_u_values)
        D_u = D_u_values(i);
        
        m = 0; % Symmetry parameter

        sol = pdepe(m, @(x,t,u,DuDx)pdeFunction(x, t, u, DuDx, D_u), @initialConditions, @boundaryConditions, x, t);

        % Extract solutions
        u_1 = sol(:,:,1);
        v_1 = sol(:,:,2);

        % Get the solution at t = maxt
        u_final = u_1(end, :);
        v_final = v_1(end, :);
        
        dx = x(2) - x(1);

        % Compute Left Riemann sums (excluding the last element)
        total_u_final_left = sum(u_final(1:end-1)) * dx;
        total_v_final_left = sum(v_final(1:end-1)) * dx;

        % Compute Right Riemann sums (excluding the first element)
        total_u_final_right = sum(u_final(2:end)) * dx;
        total_v_final_right = sum(v_final(2:end)) * dx;

        % Compute Average Riemann sums
        total_u_final_avg = (total_u_final_left + total_u_final_right) / 2;
        total_v_final_avg = (total_v_final_left + total_v_final_right) / 2;

        % Accumulate average values over seeds
        average_u_final_over_seeds(i) = average_u_final_over_seeds(i) + total_u_final_avg / L;
        average_v_final_over_seeds(i) = average_v_final_over_seeds(i) + total_v_final_avg / L;
    end
end

% Compute the average u and v values across all seeds
average_u_final_over_seeds = average_u_final_over_seeds / num_seeds;
average_v_final_over_seeds = average_v_final_over_seeds / num_seeds;

csvwrite('NP_average_uv_50.csv', [D_u_values' average_u_final_over_seeds' average_v_final_over_seeds']);

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
    
    % Compute pi_L and pi_H
    if u1 + u2 == 0
        frac = 0;
    else
        frac = u1 / (u1 + u2);
    end
    
    pi_u = (V-C)/2 * frac + V * (1 - frac);
    pi_v = 0 * frac + V/2 * (1 - frac);
    
    % Define the PDEs
    c = [1; 1];
    
    f = [D_u * Du1Dx;
         D_v * Du2Dx];
    
    s = [u1 * (pi_u - kappa * (u1 + u2));
         u2 * (pi_v - kappa * (u1 + u2))];
end



function [pl, ql, pr, qr] = boundaryConditions(xl, ul, xr, ur, t)
    % Left boundary conditions (zero flux)
    pl = [0; 0];
    ql = [1; 1];
    
    % Right boundary conditions (zero flux)
    pr = [0; 0];
    qr = [1; 1];
end


function u0 = initialConditions(x)
    u0 = [4000/9 + 4 * randn - 2; 2000/9 + 2 * randn - 1]; % Initial conditions for u, v, and n
end