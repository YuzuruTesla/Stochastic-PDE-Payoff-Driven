T = readtable('P_average_pi.csv');

w_v = T{:,1};
p_H   = T{:,2};
p_D   = T{:,3};

figure;
hold on;
dotSize = 80;
h_scatter_u = scatter(w_v, p_H, dotSize, 'filled', 'DisplayName', 'Average payoff $p_H$','MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
dotSize = 50;
h_scatter_v = scatter(w_v, p_D, dotSize, 'filled', 'DisplayName', 'Average payoff $p_D$','MarkerEdgeColor',[255/255 95/255 5/255],'MarkerFaceColor',[255/255 95/255 5/255]);
yline(2/3,  'LineWidth', 3, 'DisplayName', 'Uniform $p_H$ \& $p_D$', 'Color', [0 0 0], 'LineStyle', '--');
h_line_critical = xline(0.849, '--','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $w_v^*$');
xlim([0.75, 1.1]);
xlabel('Payoff Sensitivity of Doves $w_v$','FontSize',23,'FontName','Arial','Interpreter','latex');
ylabel('Payoff', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 20);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
ax.TickLabelInterpreter = 'latex';
grid on;
hold off;