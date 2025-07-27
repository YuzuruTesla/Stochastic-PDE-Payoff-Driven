data = load('NP_average_uv_50.csv'); 

D_u = data(:, 1);    % First column: Diffusivity of Hawks
u = data(:, 2);      % Second column: Average Density of u
v = data(:, 3);      % Third column: Average Density of v

figure;
hold on;
dotSize = 80;
h_scatter_u = scatter(D_u, u, dotSize, 'filled', 'DisplayName', 'Hawks $u$','MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
dotSize = 50;
h_scatter_v = scatter(D_u, v, dotSize, 'filled', 'DisplayName', 'Doves $v$','MarkerEdgeColor',[255/255 95/255 5/255],'MarkerFaceColor',[255/255 95/255 5/255]);
yline(4000/9,'LineWidth', 3, 'DisplayName', 'Uniform $u$', 'Color', [1 0.8 0], 'LineStyle', '--');
yline(2000/9, 'LineWidth', 3, 'DisplayName', 'Uniform $v$', 'Color', [0 0 0], 'LineStyle', '-.');
h_line_critical = xline(4.917, '--','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $D_u^*$');
ylim([100, 500]);
xlabel('Diffusivity of Hawks $D_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Density of Individuals', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 20);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
ax.TickLabelInterpreter = 'latex';
grid on;
hold off;