function [bestx, bestf, bestc, archive, ul_count] = ea_dualpopsolver(funh_obj1, funh_obj2, num_xvar, lb, ub, initmatrix, funh_con1, funh_con2, param, varargin)
% This ea solver mains two population
% create argument parser
p = inputParser;
addRequired(p, 'funh_obj1');
addRequired(p, 'funh_obj2');
addRequired(p, 'num_xvar');
addRequired(p, 'lb');
addRequired(p, 'ub');
addRequired(p, 'initmatrix');
addRequired(p, 'funh_con1');
addRequired(p, 'funh_con2');
addRequired(p, 'param');
addParameter(p, 'visualize', false);
addParameter(p, 'pf', []);
addParameter(p, 'upf', []);
addParameter(p, 'xu', []);
parse(p, funh_obj1, funh_obj2, num_xvar, lb, ub, initmatrix, funh_con1, funh_con2, param,  varargin{:});
%------- interpret argument ---------------
funh_obj1 = p.Results.funh_obj1;
funh_obj2 = p.Results.funh_obj2;
num_xvar = p.Results.num_xvar;
lb = p.Results.lb;
ub = p.Results.ub;
initmatrix= p.Results.initmatrix;
funh_con1 = p.Results.funh_con1;
funh_con2 = p.Results.funh_con2;
param = p.Results.param;
visualize = p.Results.visualize;
pf = p.Results.pf;
upf = p.Results.upf;
xu = p.Results.xu;
%-----------

record = false;
if visualize && record
    f1 = figure(2);
    obj = VideoWriter('moving.avi');
    obj.Quality= 100;
    obj.FrameRate = 25;
    open(obj);
end
if visualize
    f1 = figure(1);
end

% Initialization
% make sure initmatrix is unique
initmatrix = unique(initmatrix,'rows', 'stable');
n_init = size(initmatrix, 1);
n_rest = param.popsize - n_init;
X_pop =  unifrnd(repmat(lb,n_rest,1),repmat(ub,n_rest,1));
X_pop  = [X_pop; initmatrix];
% make sure initialization unique
X_pop = unique(X_pop, 'rows','stable');


% 2 active populations
population_ll = solutions();
population_ul = solutions();

% initialization using only LL objectives.
F1 = funh_obj1(X_pop);
C1 = funh_con1(X_pop);

population_ll.add(X_pop, [], F1, [], C1, []);
population_ul.add(X_pop, [], F1, [], C1, []);

% archive x for checking
archive_x = [];
archive_x = [archive_x; X_pop];
archive_ulx = [];

if visualize
    ulf_ulpop = funh_obj2(population_ul.xus);
    % ulc_ulpop = funh_con2(population_ul.xus);
    ulf_llpop = funh_obj2(population_ll.xus);
    % ulc_llpop = funh_con2(population_ll.xus);
    plotMO_dual(f1, population_ul.FUs, population_ll.FUs, pf, 0, upf, ulf_ulpop, ulf_llpop);
end

first_encoutner_flag = true;

for ii = 1:param.gen-1
    child_X = generate_child_DE(lb, ub, population_ul.XUs, param);
    [child_X, ~, ~] = unique(child_X, 'rows', 'stable');
    [child_X, ~] = remove_repeated_solution(archive_x, child_X);
    archive_x = [archive_x; child_X];

    % child evaluation
    F_child = funh_obj1(child_X);
    C_child = funh_con1(child_X);

    % combine 3 populations
    population_3n = solutions();
    population_3n.add(child_X, [], F_child, [], C_child, []);
    population_3n.merge(population_ll);
    population_3n.merge(population_ul);

    % temporary population for UL evaluation
    tmp_population_3n = solutions();
    tmp_population_3n.copy(population_3n);

    % evaluation to fix objectives being possible different for two
    % population, so first evaluate on ll
    tmp_3nx = population_3n.XUs;
    tmp_3nf = funh_obj1(tmp_3nx);
    tmp_3nc = funh_con1(tmp_3nx);
    population_3n.FU = {tmp_3nf};
    population_3n.FC = {tmp_3nc};

    % ND sort:
    num_nd = population_3n.nd_sort_no_reduction();
    selection_id = 1:param.popsize;
    population_ll.clear_data();

    newpop_X = population_3n.xus;
    newpop_F = population_3n.FUs;
    newpop_C = population_3n.FCs;

    newpop_X = newpop_X(selection_id, :);
    newpop_F = newpop_F(selection_id, :);
    if ~isempty(newpop_C)
        newpop_C = newpop_C(selection_id, :);
    end
    population_ll.add(newpop_X, [], newpop_F, [], newpop_C, []);

    if num_nd >= param.popsize
        % by now population_3n has been sorted;
        FU = funh_obj2(tmp_population_3n.xus);
        FC = funh_con2(tmp_population_3n.xus);
        tmp_population_3n.FU = {FU};
        tmp_population_3n.FC = {FC};
        % redo sorting on UL objective
        tmp_population_3n.nd_sort_no_reduction();

        if first_encoutner_flag
            unique_x = unique(tmp_population_3n.XUs, 'rows','stable');
            first_encoutner_flag = false;
            archive_ulx = [archive_ulx; unique_x];
        else
            unique_x = unique(tmp_population_3n.XUs, 'rows','stable');
            lia = ismember(unique_x, archive_ulx, 'rows');
            new_uniquex = unique_x(~lia, :);
            archive_ulx = [archive_ulx; new_uniquex];
        end

        newpop_X = tmp_population_3n.xus;
        newpop_F = tmp_population_3n.FUs;
        newpop_C = tmp_population_3n.FCs;

        newpop_X = newpop_X(selection_id, :);
        newpop_F = newpop_F(selection_id, :);
        if ~isempty(newpop_C)
            newpop_C = newpop_C(selection_id, :);
        end
        if ii == 17

             fl = funh_obj1(newpop_X);
            plotMO_UL(xu, newpop_F,newpop_C,fl, ii);
        end
    end

    population_ul.clear_data();
    population_ul.add(newpop_X, [], newpop_F, [], newpop_C, []);

    if visualize
        ulf_ulpop = funh_obj2(population_ul.xus);
        % ulc_ulpop = funh_con2(population_ul.xus);
        ulf_llpop = funh_obj2(population_ll.xus);
        % ulc_llpop = funh_con2(population_ll.xus);
        llf_ulpop = funh_obj1(population_ul.xus);
        plotMO_dual(f1, llf_ulpop, population_ll.FUs, pf, ii+1, upf, ulf_ulpop, ulf_llpop);
    end

    clear population_3n
    clear tmp_population_3n
end

% correct population_ul to FL
tmpx = population_ul.XUs;
tmpf = funh_obj1(tmpx);
tmpc = funh_con1(tmpx);
population_ul.FU = {tmpf};
population_ul.FC = {tmpc};

population_ul.nd_sort();


bestx = population_ul.XUs;
bestf = population_ul.FUs;
bestc = population_ul.FCs;
archive = archive_x;
ul_count = size(archive_ulx, 1);

if visualize
close(f1);
end


end

function plotMO_dual(fighn, pop_ul, pop_ll, pf, gen, upf, ulf_ulpop, ulf_llpop)
clf(fighn);
subplot(1, 2, 1);
plot(pf(:, 1), pf(:, 2), 'k.'); hold on;
scatter(pop_ul(:, 1), pop_ul(:, 2),  30, 'r', 'filled');
scatter(pop_ll(:, 1), pop_ll(:, 2), 50, 'blue');
% xlim([0, 20]);
% ylim([0, 20]);
title(num2str(gen));
grid on;

subplot(1, 2, 2);
plot(upf(:, 1), upf(:, 2), 'k.'); hold on;
scatter(ulf_ulpop(:, 1), ulf_ulpop(:, 2),  30, 'r', 'filled');
scatter(ulf_llpop(:, 1), ulf_llpop(:, 2), 50, 'blue');
title(num2str(gen));
% xlim([0, 20]);
% ylim([0, 20]);
grid on;

pause(0.2);
end


function plotMO_UL( xu, F, C, fl, gen)
prob = DS4D(3,2);
f2 = figure;
load("C:\Users\z3276872\Documents\cec24\problems\DS4m_ULPF1025.mat");

scatter(pf(:, 1), pf(:, 2), 10, [0 0.4470 0.7410], 'filled',  'DisplayName','UL Pareto front'); hold on;

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


id = C<=0;
scatter(F(id, 1), F(id, 2), 60, [0.9290 0.6940 0.1250], 'filled', 'DisplayName','LL solutions on UL objective');
legend('FontSize', 18, 'Location','northoutside');
xlabel('F_{1}', 'FontSize', 18);
ylabel('F_{2}', 'FontSize', 18,'Rotation',0);

% ax = gca;
% ax.XAxis.FontSize = 18;
% ax.YAxis.FontSize = 18;
xlim([-10, 10]);
ylim([0, 1]);


savename = sprintf("deceptiveUL.fig");
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);

savename = sprintf("deceptiveUL.png");
savename = fullfile(pwd, 'postprocess', savename);
saveas(f2, savename);
close(f2);

f2 = figure;
[FL, FC] = prob.evaluate_l(repmat(xu, num_ll, 1), ps);
scatter(FL(:, 1), FL(:, 2), 60, 'red', 'filled', 'DisplayName','LL Pareto front'); hold on;
scatter(fl(:, 1), F(:, 2), 60, [0.9290 0.6940 0.1250], 'filled', 'DisplayName','LL solutions');
box on;
grid on;

legend('FontSize', 18, 'Location','northoutside');
xlabel('f_{1}', 'FontSize', 18);
ylabel('f_{2}', 'FontSize', 18,'Rotation',0);

% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;

xlim([1, 2]);
ylim([0, 1]);


savename = sprintf("deceptiveLL.fig");
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);

savename = sprintf("deceptiveLL.png");
savename = fullfile(pwd, 'postprocess', savename);
saveas(f2, savename);
close(f2);

end