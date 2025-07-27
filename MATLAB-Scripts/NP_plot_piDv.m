data = load('NP_average_pHD_50.csv'); 

D_u = data(:, 1);      % First column: Diffusivity of Hawks
p_H = data(:, 2);      % Second column: Average payoff of u
p_D = data(:, 3);      % Third column: Average payoff of v

figure;
hold on;
dotSize = 80;
h_scatter_u = scatter(D_u, p_H, dotSize, 'filled', 'DisplayName', 'Average payoff $p_H$','MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
dotSize = 50;
h_scatter_v = scatter(D_u, p_D, dotSize, 'filled', 'DisplayName', 'Average payoff $p_D$','MarkerEdgeColor',[255/255 95/255 5/255],'MarkerFaceColor',[255/255 95/255 5/255]);
yline(2/3,  'LineWidth', 3, 'DisplayName', 'Uniform $p_H$ \& $p_D$', 'Color', [0 0 0], 'LineStyle', '--');
h_line_critical = xline(4.917, '--','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $D_u^*$');
ylim([0.6, 0.9]);
xlabel('Diffusivity of Hawks $D_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Payoff', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 20);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
hold off;