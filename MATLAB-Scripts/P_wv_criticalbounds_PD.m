function P_wv_criticalbounds_PD
    % 1) Parameters & equilibrium
    clc; close all;

    % parameters
    V     = 4;
    C     = 6;
    kappa = 0.001;

    % fully payoff-driven motion
    D_u   = 0.1;
    D_v   = 0.1;

    % % mix of diffusive and payoff effects
    % D_u   = 4.2;
    % D_v   = 0.1;

    % equilibrium from ODE
    u0 = -(V^2*(V - C)) / (2*kappa*C^2);
    v0 =  ( V *(V - C)^2 )/(2*kappa*C^2);
    fprintf('u0: %f\n', u0);
    fprintf('v0: %f\n', v0);

    % ODE Jacobian entries at (u0,v0)
    [a1,b1,a2,b2] = CVHDJacobian(V,C);
    
    % 2) sweep w_u
    w_u_vec   = linspace(0,0.5,51);
    w_v_I     = nan(size(w_u_vec));
    w_v_II    = nan(size(w_u_vec));
    w_v_IIIa  = nan(size(w_u_vec));   % <-- new bound

    for i = 1:numel(w_u_vec)
        w_u = w_u_vec(i);
    
        % (I) trace bound
        w_v_I(i) = (u0+v0)^2 ...
                   / (D_v*V*u0*v0) ...
                   * ( D_u + D_v ...
                       + D_u*w_u*((V+C)*v0*u0)/(u0+v0)^2 );
    
        % (II) a=0 bound
        w_v_II(i) = (u0+v0)^2 / (D_u*u0*v0*V) ...
                    * ( D_u + D_u*w_u*((V+C)*v0*u0)/(u0+v0)^2 );
    
        % (IIIa) new analytic bound:
        % w_v^{IIIa} = ((u0+v0)^2 / ((a1*u0 + a2*v0)*D_v*v0*V)) * ...
        %               ( a1*D_v + b2*D_u + ((b1*u0 + b2*v0)*D_u*u0*(V+C)*w_u)/(u0+v0)^2 )
        denom = (a1*u0 + a2*v0) * D_v * v0 * V;
        if abs(denom) < 1e-14
            w_v_IIIa(i) = NaN;  % avoid blow-up if denominator is ~0
        else
            inner = a1*D_v + b2*D_u ...
                    + ((b1*u0 + b2*v0) * D_u * u0 * (V + C) * w_u) / (u0 + v0)^2;
            w_v_IIIa(i) = ((u0 + v0)^2 / denom) * inner;
        end
    end
    
    % 3) Plot
    figure; hold on; grid on; box on;
    
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    plot(w_u_vec, w_v_I,     'LineWidth', 6, 'DisplayName', '$w_v^{\mathrm{I}}$',        'Color', [0 0.24 0.47]);
    plot(w_u_vec, w_v_II,    'LineWidth', 6, 'DisplayName', '$w_v^{\mathrm{II}}$',     'Color', [255/255  95/255   5/255]);
    plot(w_u_vec, w_v_IIIa,  'LineWidth', 6, 'DisplayName', '$w_v^{\mathrm{III}a}$', 'Color', [1 0.8 0], 'LineStyle', '--');
    
    xlabel('Payoff Sensitivity for Hawks $w_u$', 'FontSize', 23);
    ylabel('Payoff Sensitivity for Doves $w_v$', 'FontSize', 23);
    legend('Location', 'best', 'FontSize', 20);
    set(gca, 'FontSize', 20);
    axis square;
end

% Jacobian
function [a1,b1,a2,b2] = CVHDJacobian(V,C)
    a1 = V*(-C+V)*(C+2*V)/(2*C^2);
    b1 = -V*(V-C)^2/(C^2);
    a2 = V^3/(C^2);
    b2 = V*(V-C)*(C-2*V)/(2*C^2);
end
