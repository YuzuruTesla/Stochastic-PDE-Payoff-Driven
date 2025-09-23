function P_wv_criticalbounds_mix
    % 1) Parameters & equilibrium
    clc; close all;

    % parameters
    V     = 4;
    C     = 6;
    kappa = 0.001;

    % % fully payoff-driven motion
    % D_u   = 0.1;
    % D_v   = 0.1;

    % mix of diffusive and payoff effects
    D_u   = 4.2;
    D_v   = 0.1;

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
    w_v_IIIa  = nan(size(w_u_vec));   % new analytic bound (IIIa)
    w_v_IIIb  = nan(size(w_u_vec));   % new analytic bound (IIIb)
    w_v_III   = nan(size(w_u_vec));   % max(IIIa, IIIb)

    % handy constants
    den = (u0+v0)^2;
    K1  = u0*v0*(V + C)/den;
    K2  = u0*v0*V/den;
    K3  = u0^2*(V + C)/den;
    K4  = v0^2*V/den; %#ok<NASGU> % kept for clarity; cross-term cancels

    % coefficients for alpha = A w_v + B w_u + D
    %            and beta  = E w_v + F w_u + G
    % gamma = a1*b2 - b1*a2
    A = - D_u*D_v * K2;
    B =   D_u*D_v * K1;
    D_const = D_u*D_v;

    E =  D_v * V/den * (a1*u0*v0 + a2*v0^2);
    F = -D_u * (V + C)/den * (b2*u0*v0 + b1*u0^2);
    G = -a1*D_v - b2*D_u;

    gamma = a1*b2 - b1*a2;

    % sweep
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

        % (IIIa)
        denom_IIIa = (a1*u0 + a2*v0) * D_v * v0 * V;
        if abs(denom_IIIa) < 1e-14
            w_v_IIIa(i) = NaN;
        else
            inner = a1*D_v + b2*D_u ...
                  + ((b1*u0 + b2*v0) * D_u * u0 * (V + C) * w_u) / (u0 + v0)^2;
            w_v_IIIa(i) = ((u0 + v0)^2 / denom_IIIa) * inner;
        end

        % (IIIb) beta^2 - 4 alpha gamma > 0
        % alpha = A w_v + B w_u + D, beta = E w_v + F w_u + G
        % Solve: E^2 w_v^2 + [2E(Fw_u+G) - 4A gamma] w_v + [(Fw_u+G)^2 - 4 gamma (B w_u + D)] > 0
        E2   = E^2;
        Hw   = F*w_u + G;
        quad_b = 2*E*Hw - 4*A*gamma;
        quad_c = Hw^2 - 4*gamma*(B*w_u + D_const);

        % We take the larger root (lower bound for w_v) when E2 > 0.
        if E2 > 1e-14
            disc = quad_b^2 - 4*E2*quad_c;
            if disc >= 0
                w_v_IIIb(i) = ( (4*A*gamma - 2*E*Hw)/(2*E2) ) + sqrt(disc)/(2*E2);
            else
                % No real roots -> quadratic always positive (since E2>0)
                % => inequality holds for all w_v, so no lower bound needed.
                w_v_IIIb(i) = -Inf;  % makes max(IIIa,IIIb) ignore IIIb
            end
        else
            % Degenerate (E≈0): inequality reduces to linear in w_v:
            % (-4A gamma) w_v + quad_c > 0
            L = -4*A*gamma;
            if abs(L) < 1e-14
                % Completely w_v-independent; either always true or false.
                w_v_IIIb(i) = -Inf; % treat as no extra bound
            else
                thr = -quad_c / L;
                if L > 0
                    w_v_IIIb(i) = thr;  % lower bound: w_v > thr
                else
                    % upper bound in w_v (not useful for ">" style); ignore
                    w_v_IIIb(i) = -Inf;
                end
            end
        end
    end

    % (III) take the elementwise max of IIIa and IIIb, ignoring -Inf/NaN safely
    w_v_III = w_v_IIIa;
    use_b   = (isnan(w_v_IIIa) & ~isnan(w_v_IIIb)) | (w_v_IIIb > w_v_IIIa);
    w_v_III(use_b) = w_v_IIIb(use_b);

    % 3) Plot
    figure; hold on; grid on; box on;
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    plot(w_u_vec, w_v_I,     'LineWidth', 5, 'DisplayName', '$w_v^{\mathrm{I}}$',        'Color', [0 0.24 0.47] );
    plot(w_u_vec, w_v_II,    'LineWidth', 5, 'DisplayName', '$w_v^{\mathrm{II}}$',       'Color', [255/255  95/255   5/255]);
    % plot(w_u_vec, w_v_IIIa,  'LineWidth', 5, 'DisplayName', '$w_v^{\mathrm{III}a}$', 'Color', [0.4 0.7 1], 'LineStyle','-');
    % plot(w_u_vec, w_v_IIIb,  'LineWidth', 5, 'DisplayName', '$w_v^{\mathrm{III}b}$', 'Color', [1 0.8 0], 'LineStyle','--');
    plot(w_u_vec, w_v_III,   'LineWidth', 5, 'DisplayName', '$w_v^{\mathrm{III}}$', 'Color', [0.2 0.7 0.2],'LineStyle', '--');

    xlabel('Payoff Sensitivity for Hawks $w_u$', 'FontSize', 23);
    ylabel('Payoff Sensitivity for Doves $w_v$', 'FontSize', 23);
    % xlim([0,0.1]);
    % ylim([0,3]);
    legend('Location', 'best', 'FontSize', 18);
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
