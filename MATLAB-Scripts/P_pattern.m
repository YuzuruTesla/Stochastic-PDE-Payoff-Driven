clc; clear; close all;

seed = 10;
rng(seed,"twister");

L = 40;
maxt = 2000000;

x = linspace(0, L, 4001); 
t = linspace(0, maxt, 201); 

m = 0; % Symmetry parameter

sol = pdepe(m, @pdeFunction, @initialConditions, @boundaryConditions, x, t);

% Extract solutions
u = sol(:,:,1);
v = sol(:,:,2);

% Get the solution at t = maxt
u_final = u(end, :);
v_final = v(end, :);

% Get the solution at t = 0
u_initial = u(1, :);
v_initial = v(1, :);

% Parameters
V = 4;
C = 6;

% Compute payoffs p_H and p_D at t = maxt
frac_final = u_final ./ (u_final + v_final);
frac_final(isnan(frac_final)) = 0; % Handle division by zero

p_H_final = (V-C)/2 .* frac_final + V .* (1 - frac_final);
p_D_final =      0 .* frac_final + V/2 .* (1 - frac_final);

% Compute the left and right Riemann sums of payoffs at t=maxt
dx = x(2) - x(1);

% Left Riemann sums (excluding the last element)
total_payoff_p_H_final_left = sum(p_H_final(1:end-1) .* u_final(1:end-1)) * dx;
total_payoff_p_D_final_left = sum(p_D_final(1:end-1) .* v_final(1:end-1)) * dx;
total_u_final_left           = sum(u_final(1:end-1))                      * dx;
total_v_final_left           = sum(v_final(1:end-1))                      * dx;

% Right Riemann sums (excluding the first element)
total_payoff_p_H_final_right = sum(p_H_final(2:end) .* u_final(2:end)) * dx;
total_payoff_p_D_final_right = sum(p_D_final(2:end) .* v_final(2:end)) * dx;
total_u_final_right          = sum(u_final(2:end))                   * dx;
total_v_final_right          = sum(v_final(2:end))                   * dx;

% Average Riemann sums
total_payoff_p_H_final_avg = (total_payoff_p_H_final_left + total_payoff_p_H_final_right) / 2;
total_payoff_p_D_final_avg = (total_payoff_p_D_final_left + total_payoff_p_D_final_right) / 2;
total_u_final_avg          = (total_u_final_left           + total_u_final_right)           / 2;
total_v_final_avg          = (total_v_final_left           + total_v_final_right)           / 2;

average_payoff_p_H_final = total_payoff_p_H_final_avg / total_u_final_avg;
average_payoff_p_D_final = total_payoff_p_D_final_avg / total_v_final_avg;

fprintf('Average Riemann sum of p_H at t=maxt: %f\n', average_payoff_p_H_final);
fprintf('Average Riemann sum of p_D at t=maxt: %f\n', average_payoff_p_D_final);

% Compute the left and right Riemann sums of u and v at t=maxt
total_u_final_left  = sum(u_final(1:end-1)) * dx;
total_v_final_left  = sum(v_final(1:end-1)) * dx;
total_u_final_right = sum(u_final(2:end))   * dx;
total_v_final_right = sum(v_final(2:end))   * dx;

% Average Riemann sums of densities
total_u_final_avg = (total_u_final_left + total_u_final_right) / 2;
total_v_final_avg = (total_v_final_left + total_v_final_right) / 2;

average_u_final = total_u_final_avg / L;
average_v_final = total_v_final_avg / L;

fprintf('Average Riemann sum of u at t=maxt: %f\n', average_u_final);
fprintf('Average Riemann sum of v at t=maxt: %f\n', average_v_final);

% Plot p_H and p_D 
figure;
plot(x, p_H_final, 'LineWidth', 4, 'DisplayName', 'Patterned $p_H$', 'Color',[0 0.24 0.47]);
hold on;
plot(x, p_D_final, 'LineWidth', 4, 'DisplayName', 'Patterned $p_D$', 'Color',[255/255 95/255 5/255]);
yline(2/3,  'LineWidth', 4, 'DisplayName', 'Uniform $p_H$ \& $p_D$', 'Color',[0.5 0.5 0.5], 'LineStyle', '--');
xlabel('Position $x$', 'FontSize', 23, 'FontName', 'arial', 'Interpreter', 'latex');
ylabel('Payoff',            'FontSize', 23, 'FontName', 'arial', 'Interpreter', 'latex');
ylim([0.55 0.85]);
legend('show', 'Location', 'best', 'FontSize', 17, 'Interpreter', 'latex', 'FontName', 'arial');
set(gca, 'FontSize', 23);
grid on;

% Plot u and v
figure;
plot(x, u_final, 'LineWidth', 4, 'DisplayName', 'Patterned $u$', 'Color',[0 0.24 0.47]);
hold on;
plot(x, v_final, 'LineWidth', 4, 'DisplayName', 'Patterned $v$', 'Color',[255/255 95/255 5/255]);
yline(4000/9,'LineWidth', 4, 'DisplayName', 'Uniform $u$', 'Color',[0.5 0.5 0.5], 'LineStyle', '--');
yline(2000/9,'LineWidth', 4, 'DisplayName', 'Uniform $v$', 'Color',[0.5 0.5 0.5], 'LineStyle', '-.');
xlabel('Position $x$',           'FontSize', 23, 'FontName', 'arial', 'Interpreter', 'latex');
ylabel('Density of Individuals','FontSize', 23, 'FontName', 'arial', 'Interpreter', 'latex');
ylim([100 500]);
legend('show', 'Location', 'best', 'FontSize', 17, 'Interpreter', 'latex', 'FontName', 'arial');
set(gca, 'FontSize', 23);
grid on;

% -------------------------------------------------------------------------
function [c,f,s] = pdeFunction(~,~,U,DuDx)

    u  = U(1);   v  = U(2);
    ux = DuDx(1); vx = DuDx(2);

    % parameters
    D_u   = 4.2;
    D_v   = 0.1;
    w_u   = 0.05;
    w_v   = 0.85;
    kappa = 0.001;
    V = 4;
    C = 6;

    % payoffs pH (hawk) and pD (dove)
    denom = u+v;
    if denom < 1e-14
        frac = 0;
    else
        frac = u/denom;
    end
    pH = ((V-C)/2)*frac + V*(1-frac);
    pD = (V/2)*(1-frac);

    % their partials
    [dpH_du, dpH_dv] = partial_pi_u(u,v,V,C);
    [dpD_du, dpD_dv] = partial_pi_v(u,v,V,C);

    % spatial derivative of payoffs
    pH_x = dpH_du*ux + dpH_dv*vx;
    pD_x = dpD_du*ux + dpD_dv*vx;

    % exponential rule flux terms
    flux_u = D_u*ux - 2*D_u*w_u * ( u * pH_x );
    flux_v = D_v*vx - 2*D_v*w_v * ( v * pD_x );

    f = [flux_u; flux_v];

    s = [u*( pH - kappa*(u+v) );
         v*( pD - kappa*(u+v) )];

    % time‐derivative coefficients
    c = [1; 1];
end

function [pl, ql, pr, qr] = boundaryConditions(~,~,~,~,~)
    % Zero-flux (Neumann) boundary conditions
    pl = [0; 0]; ql = [1; 1];
    pr = [0; 0]; qr = [1; 1];
end

function U0 = initialConditions(~)
    U0 = [4000/9 + 0.4*randn - 0.2;
          2000/9 + 0.4*randn - 0.2];
end

function [dpdu, dpdv] = partial_pi_u(u,v,V,C)
    denom = u+v;
    if denom < 1e-14
        dpdu = 0; dpdv = 0; return;
    end
    N = ((V-C)/2)*u + V*v;
    dNdu = (V-C)/2; dNdv = V;
    dpdu = (dNdu*denom - N) / (denom^2);
    dpdv = (dNdv*denom - N) / (denom^2);
end

function [dpdu, dpdv] = partial_pi_v(u,v,V,C)
    denom = u+v;
    if denom < 1e-14
        dpdu = 0; dpdv = 0; return;
    end
    N = (V/2)*v;
    dNdu = 0; dNdv = V/2;
    dpdu = (dNdu*denom - N) / (denom^2);
    dpdv = (dNdv*denom - N) / (denom^2);
end
