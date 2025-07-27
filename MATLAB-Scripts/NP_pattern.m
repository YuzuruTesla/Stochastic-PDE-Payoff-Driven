clc; clear; close all;

seed = 3;
rng(seed,"twister");

L = 40;
maxt = 4e6;

x = linspace(0, L, 4001); 
t = linspace(0, maxt, 2001); 

m = 0; % Symmetry parameter

sol = pdepe(m, @pdeFunction, @initialConditions, @boundaryConditions, x, t);

% Extract solutions
u = sol(:,:,1);
v = sol(:,:,2);

% Get the solution at t = maxt
u_final = u(end, :);
v_final = v(end, :);

% Parameters
V = 4;
C = 6;

% Compute payoffs p_H and p_D at t = maxt
frac_final = u_final ./ (u_final + v_final);
frac_final(isnan(frac_final)) = 0; % Handle division by zero

p_H_final = (V-C)/2 .* frac_final + V .* (1 - frac_final);
p_D_final =      0 .* frac_final + V/2 .* (1 - frac_final);

% Grid spacing
dx = x(2) - x(1);

% Left Riemann sums of payoffs (excluding the last element)
total_payoff_p_H_final_left = sum(p_H_final(1:end-1) .* u_final(1:end-1)) * dx;
total_payoff_p_D_final_left = sum(p_D_final(1:end-1) .* v_final(1:end-1)) * dx;

% Right Riemann sums of payoffs (excluding the first element)
total_payoff_p_H_final_right = sum(p_H_final(2:end) .* u_final(2:end)) * dx;
total_payoff_p_D_final_right = sum(p_D_final(2:end) .* v_final(2:end)) * dx;

% Average Riemann sums of payoffs
total_payoff_p_H_final_avg = (total_payoff_p_H_final_left + total_payoff_p_H_final_right) / 2;
total_payoff_p_D_final_avg = (total_payoff_p_D_final_left + total_payoff_p_D_final_right) / 2;

% Compute average payoffs
average_payoff_p_H_final = total_payoff_p_H_final_avg / ((sum(u_final(1:end-1))*dx + sum(u_final(2:end))*dx)/2);
average_payoff_p_D_final = total_payoff_p_D_final_avg / ((sum(v_final(1:end-1))*dx + sum(v_final(2:end))*dx)/2);

fprintf('Average Riemann sum of payoffs of p_H at t=maxt: %f\n', average_payoff_p_H_final);
fprintf('Average Riemann sum of payoffs of p_D at t=maxt: %f\n', average_payoff_p_D_final);

% Compute the left and right Riemann sums of u and v at t=maxt
total_u_final_left  = sum(u_final(1:end-1)) * dx;
total_v_final_left  = sum(v_final(1:end-1)) * dx;
total_u_final_right = sum(u_final(2:end))   * dx;
total_v_final_right = sum(v_final(2:end))   * dx;

% Average densities
average_u_final = ((total_u_final_left  + total_u_final_right) / 2) / L;
average_v_final = ((total_v_final_left  + total_v_final_right) / 2) / L;

fprintf('Average Riemann sum of u at t=maxt: %f\n', average_u_final);
fprintf('Average Riemann sum of v at t=maxt: %f\n', average_v_final);

% Plot p_H and p_D 
figure;
plot(x, p_H_final, 'LineWidth', 4, 'DisplayName', 'Patterned $p_H$','Color',[0 0.24 0.47]);
hold on;
plot(x, p_D_final, 'LineWidth', 4, 'DisplayName', 'Patterned $p_D$','Color',[255/255 95/255 5/255]);
yline(2/3,  'LineWidth', 4, 'DisplayName', 'Uniform $p_H$ \& $p_D$', 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
xlabel('Position $x$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Payoff', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex'); 
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 23);
ylim([0.58,0.8]); % D_u = 4.93
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;

% Plot u and v
figure;
plot(x, u_final, 'LineWidth', 4, 'DisplayName', 'Patterned $u$','Color',[0 0.24 0.47]);
hold on;
plot(x, v_final, 'LineWidth', 4, 'DisplayName', 'Patterned $v$','Color',[255/255 95/255 5/255]);
yline(4000/9,'LineWidth', 4, 'DisplayName', 'Uniform $u$', 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
yline(2000/9, 'LineWidth', 4, 'DisplayName', 'Uniform $v$', 'Color', [0.5 0.5 0.5], 'LineStyle', '-.');
xlabel('Position $x$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Density of Individuals', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 23);
ylim([150,500]); % D_u = 4.93
% ylim([0,900]); % D_u = 10,19.55,40ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;



function [c, f, s] = pdeFunction(x, t, u, DuDx)
    % Parameters

    D_u = 4.93;
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