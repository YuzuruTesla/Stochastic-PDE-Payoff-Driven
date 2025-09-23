T = readtable('P_average_uv.csv');

w_v = T{:,1};
u   = T{:,2};
v   = T{:,3};

figure; hold on;

h_scatter_u = scatter(w_v, u,  80, 'filled', ...
    'DisplayName','Hawks $u$','MarkerEdgeColor',[0 0.24 0.47],'MarkerFaceColor',[0 0.24 0.47]);
h_scatter_v = scatter(w_v, v,  50, 'filled', ...
    'DisplayName','Doves $v$','MarkerEdgeColor',[255/255 95/255 5/255],'MarkerFaceColor',[255/255 95/255 5/255]);

yline(4000/9,'--','LineWidth',3, 'DisplayName','Uniform $u$','Color',[1 0.8 0]);
yline(2000/9,'-.','LineWidth',3, 'DisplayName','Uniform $v$','Color',[0 0 0]);
xline(0.849,'--','LineWidth',3, 'DisplayName','Critical $w_v^*$','Color',[0.5 0.5 0.5]);

xlim([0.75, 1.1]);
ylim([100, 500]);
xlabel('Payoff Sensitivity of Doves $w_v$','FontSize',23,'FontName','Arial','Interpreter','latex');
ylabel('Density of Individuals','FontSize',23,'FontName','Arial','Interpreter','latex');
legend('Location','best','FontSize',17,'Interpreter','latex','FontName','Arial');
grid on;
set(gca,'FontSize',20);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
ax.TickLabelInterpreter = 'latex';
hold off;
