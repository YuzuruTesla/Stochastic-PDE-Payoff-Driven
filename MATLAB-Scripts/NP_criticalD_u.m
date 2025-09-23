clear; clc; close all;

% parameters:
V     = 4; 
C     = 6; 
kappa = 0.001;
L     = 40;
Dv    = 0.1;  

%  a1, b1, a2, b2
a1 =  V * ( -C + V ) * ( C + 2*V ) / ( 2*C^2 );  
a2 = -V * ( V - C )^2 / ( C^2 );                 
b1 =  V^3 / ( C^2 );                             
b2 =  V * (V - C ) * ( C - 2*V ) / ( 2*C^2 );     

% m values
mvals   = (1:19)';
DuCrit  = nan(size(mvals));

% Compute critical Du for each integer m
for i = 1:numel(mvals)
    m  = mvals(i);
    k  = m*pi/L;
    numer  = a1*Dv*k^2 + b1*a2 - a1*b2;
    denom  = Dv*k^4 - b2*k^2;
    DuCrit(i) = numer/denom;
end

% Analytical 'm' where Dv*(mπ/L)^4 - b2*(mπ/L)^2 = 0
if b2/Dv > 0
    mcrit = L/pi * sqrt(b2/Dv);
else
    mcrit = NaN;   % no real root
    warning('b2/Dv ≤ 0, critical m is not real.');
end

figure;
scatter(mvals, DuCrit, 100, 'filled', 'MarkerFaceColor','[0 .24 .47]');
hold on;

% vertical dashed line at m^*
if ~isnan(mcrit)
    xline(mcrit, '--', ...
          'LineWidth',3, ...
          'Color',[0 0 0], ...
          'Label',sprintf('$m_c=%.3f$',mcrit), ...
          'LabelVerticalAlignment','top', ...
          'Interpreter','latex',...
          'FontSize',18);
end

grid on;
xlabel('Wavenumber $m$','Interpreter','latex', 'FontSize',22);
ylabel('Critical Diffusivity of Hawks $D_u^*(m)$','Interpreter','latex', 'FontSize',22);
ylim([0 400]);
% xlim([10 14]);
% ylim([4.5 5.5]);
set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',22, ...
    'LineWidth',1, ...
    'Box','off');
hold off;
