function nested_EA_forBLOP(seed, prob_str, checking_causal_relation)
tic;
% % create argument parser
% p = inputParser;
% addParameter(p, 'prob_str', 'DS5(1, 9)');
% addParameter(p, 'strategy', 'net4XLNDinit');
% addParameter(p, 'seed', 1);
% % parse input parameter
% prob_str = p.Results.prob_str;
% strategy = p.Results.strategy;
% seed = p.Results.seed;
%-------------------------------


% causal relation determination
% seed = 1;
rng(seed,'twister');
% checking_causal_relation = true;
prob = eval(prob_str);

[r, xl_causal_fl] = LL_variable_test(prob);

pf_file = fullfile(pwd, 'problems', strcat(prob.name, '_ULPF1025.mat'));
load(pf_file);


% load parameter
parameter_setting();
load('parameter.mat');

% archive xu for eliminating repeated solutions
archive_xu = [];

% set up output records
records = cell(1, 3);
records{1} = {}; % active population
records{2} = {}; % nd(accumulated) / generation
records{3} = []; % igd
records{4} = []; % FE counts

% run nested EA
XU = unifrnd(repmat(prob.ul_bl, parameter.UL_popsize,1),repmat(prob.ul_bu,parameter.UL_popsize,1));
population = solutions();

% process 1st population
pop_ULcount = 0;
pop_LLcount = 0;
for ii = 1: parameter.UL_popsize
    archive_xu = [archive_xu; XU(ii, :)];
    [population, single_xu_ULcount, single_xu_LLcount] = Determine_LLsolutions_forUL(XU(ii, :), population, parameter, prob, r, xl_causal_fl, checking_causal_relation);
    pop_ULcount = pop_ULcount + single_xu_ULcount;
    pop_LLcount = pop_LLcount + single_xu_LLcount;
end

tmp_population = solutions();
tmp_population.copy(population);
records{1} = [records{1}, tmp_population]; % active population
clear tmp_population;

tmp_population2 = solutions();
tmp_population2.copy(population);
tmp_population2.nd_sort();  % No keeping all XL for ND front xu
records{2} = [records{2}, tmp_population2];  % nd (accumulated) / generation
clear tmp_population2

current_nd_solutions = records{2}(end);
ndFU = current_nd_solutions.FUs;
igd = mean(min(pdist2(pf, ndFU),[],2));
records{3} = [records{3}, igd];

records{4} = [records{4}; pop_ULcount, pop_LLcount];


visualizationND = true;
if visualizationND
    scatter(pf(:, 1), pf(:, 2), 20, 'k', 'filled'); hold on;
    grid on;
    box on;
    scatter(ndFU(:, 1), ndFU(:, 2), 50, 'red', 'filled');
end


for jj = 1: parameter.UL_gensize

    % Generate child population
    param_tmp.popsize = parameter.UL_popsize;
    childXU = generate_child_DE(prob.ul_bl, prob.ul_bu, population.xus, param_tmp);
    childXU = unique(childXU, 'rows', 'stable');
    childXU = remove_repeated_solution(archive_xu, childXU);
    archive_xu = [archive_xu; childXU];

    % evaluate child population
    pop_ULcount = 0;
    pop_LLcount = 0;
    for ii = 1:size(childXU, 1)
        [population, single_xu_ULcount, single_xu_LLcount] = Determine_LLsolutions_forUL(childXU(ii, :), population, parameter, prob, r, xl_causal_fl, checking_causal_relation);
        pop_ULcount = pop_ULcount + single_xu_ULcount;
        pop_LLcount = pop_LLcount + single_xu_LLcount;
    end

    % save ND front, has to be done here
    tmp_population2 = solutions();
    tmp_population2.copy(records{2}(end));
    tmp_2Npopulation = solutions();
    tmp_2Npopulation.copy(population);
    tmp_population2.merge(tmp_2Npopulation);
    tmp_population2.nd_sort();                   % No keeping all XL for ND front xu
    records{2} = [records{2}, tmp_population2];  % nd (accumulated) / generation

    % ND sort
    population.DSS_newpopulation(parameter.UL_popsize, prob);

    % save active population 
    tmp_population = solutions();
    tmp_population.copy(population);
    records{1} = [records{1}, tmp_population]; % active population    
    clear tmp_population;
    clear tmp_population2
 
    % save igd
    current_nd_solutions = records{2}(end);
    ndFU = current_nd_solutions.FUs;
    igd = mean(min(pdist2(pf,ndFU),[],2));
    records{3} = [records{3}, igd];

    if visualizationND
        clf;
        scatter(pf(:, 1), pf(:, 2), 20, 'k', 'filled'); hold on;
        scatter(ndFU(:, 1), ndFU(:, 2), 50, 'red', 'filled');
        grid on;
        box on;
        t = sprintf("generation %d ", jj + 1);
        title(t)
        pause(1);
    end

    records{4} = [records{4}; pop_ULcount, pop_LLcount];
end

% record finally ND front
result_folder = fullfile(pwd, "results");
if ~exist(result_folder, "dir")
    mkdir(result_folder);
end

problem_folder = fullfile(result_folder, prob.name);
if ~exist(problem_folder, "dir")
    mkdir(problem_folder);
end

if checking_causal_relation
    save_name = sprintf("%s_CRchecking_seed_%d.mat", prob.name, seed);
else
    save_name = sprintf("%s_EA_seed_%d.mat", prob.name, seed);
end
save_file = fullfile(pwd, 'results', prob.name, save_name);
save(save_file, 'records');
toc;
end



function [r, xl_causal_fl] = LL_variable_test(prob)
% causal relation determination
r = false;

xu = unifrnd(prob.ul_bl, prob.ul_bu);
xl = unifrnd(prob.ll_bl, prob.ll_bu);


nl = length(xl);
delta = 0.1;

base_fl = prob.evaluate_l(xu, xl);
xl_causal_fl = ones(1, nl);

% test LL variables on fl
for ii = 1:nl
    delta_xl = xl(ii) + delta;

    if delta_xl >= prob.ll_bu(ii)
        delta_xl = xl(ii) - delta;
    end

    tmp_xl = xl;
    tmp_xl(ii) = delta_xl;

    tmp_fl = prob.evaluate_l(xu, tmp_xl);

    if abs(tmp_fl - base_fl) < 1e-7 % no causal relation
        xl_causal_fl(ii) = 0;
    end
end

if sum(xl_causal_fl) < nl
    r = true;
end

end




