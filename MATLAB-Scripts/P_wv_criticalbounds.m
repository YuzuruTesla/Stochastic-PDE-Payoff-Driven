function P_wv_criticalbounds
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

    L     = 40;   
    mMax  = 1000000;

    % equilibrium from ODE
    u0 = -(V^2*(V - C)) / (2*kappa*C^2);
    v0 =  ( V *(V - C)^2 )/(2*kappa*C^2);
    fprintf('u0: %f\n', u0);
    fprintf('v0: %f\n', v0);

    % ODE Jacobian entries at (u0,v0)
    [a1,b1,a2,b2] = CVHDJacobian(V,C);
    
    % 2) sweep w_u
    w_u_vec = linspace(0,0.5,51);
    w_v_I   = nan(size(w_u_vec));
    w_v_II  = nan(size(w_u_vec));
    w_v_III = nan(size(w_u_vec));
    % w_v_IV  = nan(size(w_u_vec));
    
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
    
        % (III) det=0 bound
        sols = [];
        for m = 1:mMax
            k2 = (m*pi/L)^2;
            wv = solveDetWv(a1,b1,a2,b2,u0,v0,V,C,...
                            D_u,D_v,w_u,k2);
            if isreal(wv) && wv>0
                sols(end+1) = wv;
            end
        end
        if ~isempty(sols)
            w_v_III(i) = min(sols);
        end
    
        % % (IV) \beta=0 bound
        % num = (a1*D_v + b2*D_u)*(u0+v0)^2 ...
        %       + D_u*w_u*(V+C)*u0*(b2*v0 + a2*u0);
        % den = D_v * V * v0 * (a1*u0 + b1*v0);
        % w_v_IV(i) = num ./ den;
    end
    
    % 3) Plot
    figure; hold on; grid on; box on;
    
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    plot(w_u_vec, w_v_I,   'LineWidth', 6, 'DisplayName', '$w_v^I$',   'Color', [255/255  95/255   5/255]);
    plot(w_u_vec, w_v_II,  'LineWidth', 6, 'DisplayName', '$w_v^{II}$','Color', [0 0.24 0.47]);
    plot(w_u_vec, w_v_III, 'LineWidth', 6, 'DisplayName', '$w_v^{III}$','Color', [1 0.8 0], 'LineStyle', '--');
   % plot(w_u_vec, w_v_IV,  'LineWidth', 6, 'DisplayName', '$w_v^{IV}$', 'Color', [0.2 0.7 0.2], 'LineStyle', '--');
    
    xlabel('Payoff-driven weight for Hawks $w_u$', 'FontSize', 23);
    ylabel('Payoff-driven weight for Doves $w_v$', 'FontSize', 23);
    % ylim([0,3]);  % mix of diffusive and payoff effects
    
    legend('Location', 'best', 'FontSize', 20);
    set(gca, 'FontSize', 20);
    axis square;


end

% Jacobian
function [a1,b1,a2,b2] = CVHDJacobian(V,C)
    a1 = V*(-C+V)*(C+2*V)/(2*C^2);
    a2 = -V*(V-C)^2/(C^2);
    b1 = V^3/(C^2);
    b2 = V*(V-C)*(C-2*V)/(2*C^2);
end

% solves det(A)=0 for w_v
function wv = solveDetWv(a1,b1,a2,b2,...
                         u0,v0,V,C,...
                         D_u,D_v,w_u,k2)
    A1 = a1 ...
         - D_u*k2 ...
         - D_u*w_u*(V + C)*u0*v0/(u0+v0)^2 * k2;

    B1 = b1 ...
         + D_u*w_u*(V + C)*u0^2   /(u0+v0)^2 * k2;

    A2_0    = a2;
    A2_coeff = - D_v * V * v0^2 /(u0+v0)^2 * k2;

    B2_0    = b2 - D_v*k2;
    B2_coeff =   D_v * V * u0*v0 /(u0+v0)^2 * k2;

    % det(A)= A1*(B2_0 + B2_coeff*w_v) - B1*(A2_0 + A2_coeff*w_v) = 0
    num   = -( A1*B2_0   - B1*A2_0   );
    denom =   A1*B2_coeff - B1*A2_coeff;

    if abs(denom) < 1e-14
        wv = NaN;
    else
        wv = num/denom;
    end
end
