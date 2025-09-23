V     = 4; 
C     = 6; 
kappa = 0.001;
L     = 40;
D_u   = 0.1;
D_v   = 0.1;
w_u   = 1;
w_v   = 10;
rho   = 0.01;

a_1 =  V * ( -C + V ) * ( C + 2*V ) / ( 2*C^2 );  
a_2 = -V * ( V - C )^2 / ( C^2 );                 
b_1 =  V^3 / ( C^2 );                             
b_2 =  V * (V - C) * ( C - 2*V ) / ( 2*C^2 );     

u_0 = -(V^2*(V - C)) / (2*kappa*C^2);
v_0 =  ( V *(V - C)^2 )/(2*kappa*C^2);

m_1 = -  (V+C)*v_0 / (2*(u_0+v_0)^2);
n_1 =   (V+C)*u_0 / (2*(u_0+v_0)^2);
m_2 =   V *v_0 / (2*(u_0+v_0)^2);
n_2 = -  V*u_0 / (2*(u_0+v_0)^2);

m_values = linspace(1, 4000, 2000);
k = 2 * m_values * pi / L;

% -- common term --
S = (k .* sin(k * rho)) / rho;

% -- determinant terms --
term1 = -D_u * k.^2 + a_1 + 2*D_u*w_u * u_0 * m_1 .* S;
term2 = -D_v * k.^2 + b_2 + 2*D_v*w_v * v_0 * n_2 .* S;
term3 =  b_1           + 2*D_u*w_u * u_0 * n_1 .* S;
term4 =  a_2           + 2*D_v*w_v * v_0 * m_2 .* S;
det_B = term1 .* term2 - term3 .* term4;

% -- trace --
tr_B = term1 + term2;


set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');

% -- plot det(A) --
figure;
plot(m_values, det_B, '.', 'MarkerSize', 12,'MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
yline(0, '--','Color',[1 0.8 0], 'LineWidth', 4);
grid on;
xlabel('Wavenumber $m$', 'Interpreter', 'latex');
ylabel('$\det(B)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 23);


% -- plot tr(A) --
figure;
plot(m_values, tr_B, '.', 'MarkerSize', 12,'MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
yline(0, '--','Color',[1 0.8 0], 'LineWidth', 4);
grid on;
xlabel('Wavenumber $m$', 'Interpreter','latex');
ylabel('$\mathrm{tr}(B)$', 'Interpreter','latex');
set(gca, 'FontSize', 23);
