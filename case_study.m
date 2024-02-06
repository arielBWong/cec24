prob = DS4(5, 4);

rng(1, 'twister');
xu = 1 + rand(1, 1);

% 
% f1 = figure(1);
% load("C:\Users\z3276872\Documents\cec24\problems\DS4m_ULPF1025.mat");
% 
% scatter(pf(:, 1), pf(:, 2), 10, 'k', 'filled',  'DisplayName','UL pareto front'); hold on;
% 
% num_ll = 100;
% ps = prob.PS_LL(num_ll, xu);
% [FU, FC] = prob.evaluate_u(repmat(xu, num_ll, 1), ps);
% idx = FC <= 0;
% FU_feasible = FU(idx, :);
% FU_infeasible = FU(~idx, :);
% 
% scatter(FU_feasible(:, 1), FU_feasible(:, 2), 60, 'red', 'filled', 'DisplayName','Feasible LL PS solutions');
% scatter(FU_infeasible(:, 1), FU_infeasible(:, 2), 10, 'red', 'DisplayName','Infeasible LL PS solutions');
% box on;
% grid on;
% 
% 
% box on;
% grid on;


% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% 
% xlabel('F_{1}', 'FontSize', 20);
% ylabel('F_{2}', 'FontSize', 20,'Rotation',0);
% legend('FontSize', 15, 'Location','southwest');
% 
% savename = sprintf("%s_UL_pf_and_llps", prob.name);
% savename = fullfile(pwd, 'postprocess', savename);
% savefig(savename);
% close(f1);


llpf = prob.PF_LL(513, xu);

% figure('Position',[100, 100, 600, 600]);
scatter(llpf(:, 1), llpf(:, 2), 10, 'r', 'filled', 'DisplayName', 'LL pareto front'); hold on;
box on;
grid on;




population = solutions();
parameter.UL_popsize = 50;
parameter.LL_popsize = 50;

parameter.UL_gensize = 200;
parameter.LL_gensize = 200;
[population, single_xu_ULcount, single_xu_LLcount] = Determine_LLsolutions_forUL(xu, population, parameter, prob, 0, [], 0);

% F = population.FLs;
% 
% scatter(F(:, 1), F(:, 2), 60, [0.4660 0.6740 0.1880], 'filled', 'DisplayName', 'Search results');
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;
% 
% 
% xlabel('f_{1}', 'FontSize', 20);
% ylabel('f_{2}', 'FontSize', 20,'Rotation',0);
% legend('FontSize', 15, 'Location','northeast');
% 
% savename = sprintf("%s_LL_pf_and_llsearch", prob.name);
% savename = fullfile(pwd, 'postprocess', savename);
% savefig(savename);
% close();

F = population.FUs;
scatter(F(:, 1), F(:, 2), 60, [0.4660 0.6740 0.1880], 'filled', 'DisplayName', 'LL search results on UL');
a = 0;
legend('FontSize', 15);
xlim([0, 30]);
ylim([0, 30]);
savename = sprintf("%s_UL_pf_and_llsearch_part", prob.name);
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);
close();
