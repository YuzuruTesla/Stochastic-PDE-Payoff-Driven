clear; clc; close all;

% 1) Payoff parameters
V = 4; 
C = 6;  % C>V>0

% The Jacobian entries at the interior equilibrium:
a1 =  V*( -C + V )*( C + 2*V ) / (2*C^2);
a2 = -V*( V - C )^2 / ( C^2 );
b1 =  V^3 / ( C^2 );
b2 =  V*( V - C )*( C - 2*V ) / (2*C^2 );

% 2) Create a 2D grid of (D_v, D_u)
Nv     = 201;          
Nu     = 301;          
Dv_max = 0.2;            
Du_max = 3;            

Dv_vals = linspace(0, Dv_max, Nv);
Du_vals = linspace(0, Du_max, Nu);

[Dv_grid, Du_grid] = meshgrid(Dv_vals, Du_vals);

% 3) Evaluate both conditions

% Condition 1: D_u*b2 + D_v*a1 > 0
cond1 = Du_grid.*b2 + Dv_grid.*a1;  

% Condition 2: (a1*Dv + b2*Du)^2 - 4 Du Dv (a1 b2 - a2 b1) > 0
cond2 = ( a1.*Dv_grid + b2.*Du_grid ).^2 ...
           - 4.*Du_grid.*Dv_grid.*( a1*b2 - a2*b1 );

% 4)
mask1  = (cond1>0);
mask2  = (cond2>0);
phase_map  = 2*mask1 + mask2; 

% 5) Plot
figure;
p = pcolor(Dv_grid, Du_grid, phase_map);
shading flat; 


my_colormap = [
    0    0.24   0.47;    % 0        (neither)
    255/255  95/255   5/255; % 1    (Condition 2 only)
    0.6196 0.8196 0.4824;    % 2    (Condition 1 only)
    1     0.8   0   % 3             (both)
];
colormap(my_colormap);

caxis([0 3]);

xlabel('Diffusivity of Doves $D_v$','FontSize',22,'Interpreter', 'latex');
ylabel('Diffusivity of Hawks $D_u$','FontSize',22,'Interpreter', 'latex');

hold on;

c1 = contour(Dv_grid, Du_grid, cond1, [0 0], 'LineColor','k','LineWidth',2);
c2 = contour(Dv_grid, Du_grid, cond2, [0 0], 'LineColor','k','LineWidth',2);


h_none   = plot(nan, nan, 's', ...
    'MarkerFaceColor', my_colormap(1,:), ...
    'MarkerEdgeColor', my_colormap(1,:), ...
    'MarkerSize', 40);

h_cond2  = plot(nan, nan, 's', ...
    'MarkerFaceColor', my_colormap(2,:), ...
    'MarkerEdgeColor', my_colormap(2,:), ...
    'MarkerSize', 40);

h_cond1  = plot(nan, nan, 's', ...
    'MarkerFaceColor', my_colormap(3,:), ...
    'MarkerEdgeColor', my_colormap(3,:), ...
    'MarkerSize', 40);

h_both   = plot(nan, nan, 's', ...
    'MarkerFaceColor', my_colormap(4,:), ...
    'MarkerEdgeColor', my_colormap(4,:), ...
    'MarkerSize', 40);

legendHandles = [h_none, h_cond2, h_cond1, h_both];
legendLabels = {
    'Neither condition', ...
    'Condition 2 only', ...
    'Condition 1 only', ...
    'Both conditions'
};

legend(legendHandles, legendLabels, ...
    'Interpreter','LaTeX', ... 
    'Location','bestoutside', ...
    'FontSize',20, ...
    'Box','off');

ax = gca;                    
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
ax.XTick = 0 : 0.05 : Dv_max; 
ax.YTick = 0 : 0.5  : Du_max;
hold off;



