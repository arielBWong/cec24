prob = DS4D(5, 4);

rng(2, 'twister');
xu = 1 + rand(1, 1);


f1 = figure(1);
load("C:\Users\z3276872\Documents\cec24\problems\DS4m_ULPF1025.mat");

scatter(pf(:, 1), pf(:, 2), 10, [0.4660 0.6740 0.1880], 'filled',  'DisplayName','UL Pareto front'); hold on;

num_ll = 100;
ps = prob.PS_LL(num_ll, xu);
[FU, FC] = prob.evaluate_u(repmat(xu, num_ll, 1), ps);
idx = FC <= 0;
FU_feasible = FU(idx, :);
FU_infeasible = FU(~idx, :);

scatter(FU_feasible(:, 1), FU_feasible(:, 2), 60, 'red', 'filled', 'DisplayName','Feasible LL PS solutions');
scatter(FU_infeasible(:, 1), FU_infeasible(:, 2), 10, 'red', 'DisplayName','Infeasible LL PS solutions');
box on;
grid on;


box on;
grid on;