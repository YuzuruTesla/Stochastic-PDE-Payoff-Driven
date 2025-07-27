clc; clear; close all;

L     = 40;               
maxt  = 4e6;           
x     = linspace(0, L, 4001);
t     = linspace(0, maxt, 201);
V     = 4;                % Hawk–Dove game parameters
C     = 6;
w_v_values = 0.75 : 0.005 : 1.1;


num_seeds = 50;

N = numel(w_v_values);
avg_u = zeros(1, N);
avg_v = zeros(1, N);
avg_pH = zeros(1, N);
avg_pD = zeros(1, N);

dx = x(2)-x(1);

for seed = 1:num_seeds
    rng(seed, 'twister');
    for i = 1:N
        w_v = w_v_values(i);
        sol = pdepe(0, ...
            @(x,t,U,DuDx) pdeFunction(x,t,U,DuDx,w_v), ...
            @initialConditions, @boundaryConditions, x, t);
        
        u = sol(:,:,1);  v = sol(:,:,2);
        u_f = u(end,:);  v_f = v(end,:);
        
        % Population densities
        U_left  = sum(u_f(1:end-1))   * dx;
        U_right = sum(u_f(2:end))     * dx;
        V_left  = sum(v_f(1:end-1))   * dx;
        V_right = sum(v_f(2:end))     * dx;
        
        total_U = (U_left + U_right)/2;   % integrated u
        total_V = (V_left + V_right)/2;   % integrated v
        
        avg_u(i) = avg_u(i) + total_U / L;   % density of u
        avg_v(i) = avg_v(i) + total_V / L;   % density of v
        
        % Payoffs at t=maxt
        frac = u_f ./ (u_f + v_f);
        frac(isnan(frac)) = 0;
        pH = ((V-C)/2) .* frac + V .* (1 - frac);
        pD = (V/2)     .* (1 - frac);
        
        % Weighted Riemann sums for payoffs
        PHL = sum(pH(1:end-1) .* u_f(1:end-1)) * dx;
        PHR = sum(pH(2:end)   .* u_f(2:end))   * dx;
        PDL = sum(pD(1:end-1) .* v_f(1:end-1)) * dx;
        PDR = sum(pD(2:end)   .* v_f(2:end))   * dx;
        
        total_pH = (PHL + PHR)/2;
        total_pD = (PDL + PDR)/2;
        
        avg_pH(i) = avg_pH(i) + total_pH / total_U;   % avg payoff per Hawk
        avg_pD(i) = avg_pD(i) + total_pD / total_V;   % avg payoff per Dove
    end
end

% Average over seeds
avg_u  = avg_u  ./ num_seeds;
avg_v  = avg_v  ./ num_seeds;
avg_pH = avg_pH ./ num_seeds;
avg_pD = avg_pD ./ num_seeds;

% Write two CSVs
T_uv   = table(w_v_values', avg_u', avg_v', 'VariableNames', {'w_v','u_density','v_density'});
T_pay  = table(w_v_values', avg_pH', avg_pD','VariableNames', {'w_v','avg_pH','avg_pD'});

writetable(T_uv,  'P_average_uv_test.csv');
writetable(T_pay, 'P_average_pi_test.csv');


% PDE function

function [c,f,s] = pdeFunction(~, ~, U, DuDx, w_v)

    u  = U(1);   v  = U(2);
    ux = DuDx(1); vx = DuDx(2);

    % parameters
    D_u   = 4.2;
    D_v   = 0.1;
    w_u   = 0.05;
    kappa = 0.001;
    V = 4;
    C = 6;

    % payoffs p_H (hawk) and p_D (dove)
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

    % --- exponential rule f(x)=exp(x) ⇒ f'/f = 1
    % flux terms (so that ∂_x flux gives diffusion minus movement term)
    flux_u = D_u*ux ...
           - 2*D_u*w_u * ( u * pH_x );
    flux_v = D_v*vx ...
           - 2*D_v*w_v * ( v * pD_x );

    f = [flux_u;
         flux_v];

    s1 = u*( pH - kappa*(u+v) );
    s2 = v*( pD - kappa*(u+v) );
    s  = [s1; s2];

    % time‐derivative coefficients
    c = [1; 1];
end


function [pl, ql, pr, qr] = boundaryConditions(~, Ul, ~, Ur, ~)
    % Zero-flux (Neumann) boundary conditions => flux=0 => pl=0, ql=1
    pl = [0; 0];
    ql = [1; 1];
    pr = [0; 0];
    qr = [1; 1];
end

function U0 = initialConditions(~)
    U0 = [ 4000/9 + 0.4*randn - 0.2; 
           2000/9 + 0.4*randn - 0.2 ];
end

% partial_pi_u(u,v): partial derivatives for pi_u
function [dpdu, dpdv] = partial_pi_u(u,v,V,C)
    % pi_u = [ (V-C)/2 * u + V*v ] / (u+v)
    denom = u+v;
    if denom<1e-14
        dpdu = 0; dpdv=0;
        return;
    end
    N = ((V-C)/2)*u + V*v;  % numerator
    dNdu = (V-C)/2;
    dNdv = V;

    % pi_u= N/den => partial wrt u => (dNdu*den - N* dden/du)/den^2
    % dden/du=1
    dpdu = ( dNdu*denom - N*1 ) / ( denom^2 );
    dpdv = ( dNdv*denom - N*1 ) / ( denom^2 );
end

% partial_pi_v(u,v): partial derivatives for pi_v
function [dpdu, dpdv] = partial_pi_v(u,v,V,C)
    % For your hawk-dove, pi_v= (V/2)*( v/(u+v) ) ignoring c=0
    denom = u+v;
    if denom<1e-14
        dpdu=0; dpdv=0;
        return;
    end
    N = (V/2)*v; 
    dNdu=0;  
    dNdv= V/2;

    dpdu = ( dNdu*denom - N*1 ) / ( denom^2 );
    dpdv = ( dNdv*denom - N*1 ) / ( denom^2 );
end