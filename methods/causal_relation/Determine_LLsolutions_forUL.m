function [population, UL_count, LL_count] = Determine_LLsolutions_forUL(xu, population, parameter, prob, r, xl_causal_fl, checking_causal_relation)
% xl_causal_fl: 1 - related to fl, 0 - not related to fl

visualization_ULadd = false;
visualization_UL = false;
pf = [];

% fix XL corresponding UL related variables, if necessary
if r && checking_causal_relation == 1
    % for normal search only search for partial variables
    related_to_fl = true;
    [XL, LL_count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl);
    UL_count = 0;
    if visualization_UL
        figure('Position', [100, 100, 600, 600]);
        tmp_fu = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);
        LL_ps = prob.PS_LL(129, xu);
        pf = prob.evaluate_u(repmat(xu, 129, 1), LL_ps);
        scatter(pf(:, 1), pf(:, 2), 20, 'k', 'filled'); hold on;
        scatter(tmp_fu(:, 1), tmp_fu(:, 2), 30, 'red', 'filled');
    end

    refined_XL = [];
    for jj = 1:size(XL, 1)
        xl = XL(jj, :);
        related_to_fl = false;
        [XL_update, ul_count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, xl);

        if visualization_ULadd
            LL_ps = prob.PS_LL(129, xu);
            pf = prob.evaluate_u(repmat(xu, 129, 1), LL_ps);
        end
        UL_count = UL_count + ul_count;       
        refined_XL = [refined_XL; XL_update]; 
    end

    XL = refined_XL;
    [FL, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);
    LL_count = LL_count + size(XL, 1);

elseif r && checking_causal_relation == 2
    % for normal LL search only search partial variables 
    related_to_fl = true;
    [XL, LL_count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, []);

    if visualization_UL
        % put random number in UL related values
        num_xl = size(XL, 1);
        v = xl_causal_fl == 0;
        LL_nvar = sum(v);
        sub_ub = prob.ll_bu(v);
        sub_lb = prob.ll_bl(v);
        X_pop =  unifrnd(repmat(sub_lb,num_xl,1),repmat(sub_ub,num_xl,1));
        tmpXL = XL;
        tmpv = repmat(v, num_xl, 1);
        tmpXL(tmpv) = X_pop;

        figure('Position',[50, 50, 600, 600]);
        % [tmp_fu, tmp_fc] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), tmpXL);
        [tmp_fu, tmp_fc] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);

        feasible_ids = tmp_fc <= 0 ;
        tmp_fu = tmp_fu(feasible_ids, :);
        LL_ps = prob.PS_LL(129, xu);
        pf = prob.evaluate_u(repmat(xu, 129, 1), LL_ps);
        scatter(pf(:, 1), pf(:, 2), 20, 'k', 'filled'); hold on;
        scatter(tmp_fu(:, 1), tmp_fu(:, 2), 30, 'red', 'filled');
    end
    
    % for partial UL search
    v = xl_causal_fl == 1;
    related_ids = repmat(v, size(XL, 1), 1);
    XL_LLpartial = XL(related_ids);
    XL_LLpartial = reshape(XL_LLpartial, size(XL, 1), sum(v));

    funh_obj = @(x)ULsearch_objective(xu, x, prob);
    funh_con = @(x)ULsearch_constraint(xu, x, prob);
    v = xl_causal_fl == 0;
    LL_nvar = sum(v);
    sub_ub = prob.ll_bu(v);
    sub_lb = prob.ll_bl(v);
    p.gen = parameter.LL_gensize;
    p.popsize = parameter.LL_popsize;

    if visualization_ULadd
        LL_ps = prob.PS_LL(129, xu);
        pf = prob.evaluate_u(repmat(xu, 129, 1), LL_ps);
    end

    [XL, ~, ~, archive, ~] = ea_solver(funh_obj, LL_nvar, sub_lb, sub_ub, [], funh_con, p, 'visualize', visualization_ULadd, 'pf', pf, 'termination_criterion', 0,  'init_addon', XL_LLpartial, 'xposition_indicator', v);
    UL_count = size(archive, 1);
    [FL, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);
    LL_count = LL_count + size(XL, 1);

else
    % normal LL search
    [XL, FL, FLC, LL_count] = LL_Evolution(xu, parameter, prob, [], false, 0);
    UL_count = 0;
end

if visualization_UL
    [tmp_fu, tmp_fc] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);
    feasible_ids = tmp_fc <= 0 ;
    tmp_fu = tmp_fu(feasible_ids, :);
    scatter(tmp_fu(:, 1), tmp_fu(:, 2), 50, 'blue', 'LineWidth', 2);
    box on;
    grid on;
    pause(1);
    close(f2);
end

% continue on UL evolution
[FU, FC] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);
UL_count = UL_count + size(XL, 1);

population.add(xu, XL, FU, FL, FC, FLC);

end

function [XL, count] = partial_search(xu, xl_causal_fl, prob, parameter, related_to_fl, xl)
v = xl_causal_fl == related_to_fl;
sub_ub = prob.ll_bu(v);
sub_lb = prob.ll_bl(v);

p.gen = parameter.LL_gensize;
p.popsize = parameter.LL_popsize;
if related_to_fl
    funh_obj = @(x)partial_search_objective_LL(x, xu, v, prob);
    funh_con = @(x)partial_search_constraint_LL(x, xu, v, prob);
else
    funh_obj = @(x)partial_search_objective_UL(x, xu, xl, v, prob);
    funh_con = @(x)partial_search_constraint_UL(x, xu, xl, v, prob);
end

LL_nvar = sum(v);
[XL_partial, ~, ~, archive, ~] = gsolver(funh_obj, LL_nvar, sub_lb, sub_ub, [], funh_con, p, 'visualize', false, 'pf', [], 'termination_criterion', 0);

count = size(archive.sols, 1);
if related_to_fl
    tmp_xl = zeros(size(XL_partial, 1), length(prob.ll_bu));
    tmp_v = repmat(v, size(XL_partial, 1), 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
else
    additional_num = size(XL_partial, 1);
    tmp_v = repmat(v, additional_num, 1);
    tmp_xl = repmat(xl, additional_num, 1);
    tmp_xl(tmp_v) = XL_partial;
    XL = tmp_xl;
end

end 

function [f] = partial_search_objective_LL(x, xu, v, prob)
xl = zeros(1, length(prob.ll_bu));
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[f, ~] = prob.evaluate_l(tmp_xu, tmp_xl);
end

function [c] = partial_search_constraint_LL(x, xu, v, prob)
xl = zeros(1, length(prob.ll_bu));
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[~, c] = prob.evaluate_l(tmp_xu, tmp_xl);
end


% Partial search evaluation objective function
function [f] = partial_search_objective_UL(x, xu, xl, v, prob)
% v indicate LL variables relating to UL objective
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[f, ~] = prob.evaluate_u(tmp_xu, tmp_xl);
end

function [c] = partial_search_constraint_UL(x, xu, xl, v, prob)
[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[~, c] = prob.evaluate_u(tmp_xu, tmp_xl);
end

function[f] = ULsearch_objective(xu, xl, prob)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_u(xu, xl);
end

function[c] = ULsearch_constraint(xu, xl, prob)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_u(xu, xl);
end

function[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v)
num_x = size(x, 1);
tmp_v = repmat(v, num_x, 1);

tmp_xu = repmat(xu, num_x, 1);
tmp_xl = repmat(xl, num_x, 1);
tmp_xl(tmp_v) = x;
end