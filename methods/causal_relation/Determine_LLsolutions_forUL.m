function [population, UL_count, LL_count] = Determine_LLsolutions_forUL(xu, population, parameter, prob, r, xl_causal_fl)
% xl_causal_fl: 1 - related to fl, 0 - not related to fl

visualization = false;
[XL, FL, FLC, LL_count] = LL_Evolution(xu, parameter, prob, [], false, 0);

UL_count = 0;

if visualization
    f1 = figure(1);
    [vis_fu, vis_cu] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);
    scatter(vis_fu(:, 1), vis_fu(:, 2), 50, 'r', 'filled'); hold on;
end

% fix XL corresponding UL related variables, if necessary
if r
    v = xl_causal_fl == 0;
    sub_ub = prob.ll_bu(v);
    sub_lb = prob.ll_bl(v);
    refined_XL = [];

    for jj = 1:size(XL, 1)
        xl = XL(jj, :);

        p.gen = parameter.LL_gensize;
        p.popsize = parameter.LL_popsize;

        % how to run experiments on UL? use binary v to assign value back
        funh_obj = @(x)partial_search_objective(x, xu, xl, v, prob);
        funh_con = @(x)partial_search_constraint(x, xu, xl, v, prob);

        LL_nvar = sum(v);
        [XL_partial, ~, ~, ~, ul_count] = gsolver(funh_obj, LL_nvar, sub_lb, sub_ub, [], funh_con, p, 'visualize', false, 'pf', [], 'termination_criterion', 0);
        UL_count = UL_count + ul_count;

        additional_num = size(XL_partial, 1);
        tmp_v = repmat(v, additional_num, 1);
        tmp_xl = repmat(xl, additional_num, 1);
        tmp_xl(tmp_v) = XL_partial;
        
        refined_XL = [refined_XL; tmp_xl];      
    end

    XL = refined_XL;
else
    % no action
end

% continue on UL evolution

[FU, FC] = prob.evaluate_u(repmat(xu, size(XL, 1), 1), XL);
[FC, FLC] = prob.evaluate_l(repmat(xu, size(XL, 1), 1), XL);

UL_count = UL_count + size(XL, 1);
LL_count = LL_count + size(XL, 1);

population.add(xu, XL, FU, FL, FC, FLC);

if visualization
    scatter(FU(:, 1), FU(:, 2), 60, 'b');
end
end


% Partial search evaluation objective function
function [f] = partial_search_objective(x, xu, xl, v, prob)
% v indicate LL variables relating to UL objective

[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[f, ~] = prob.evaluate_u(tmp_xu, tmp_xl);
end

function [c] = partial_search_constraint(x, xu, xl, v, prob)

[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v);
[~, c] = prob.evaluate_u(tmp_xu, tmp_xl);
end

function[tmp_xu, tmp_xl] = assign_partial_xl(x, xu, xl, v)
num_x = size(x, 1);
tmp_v = repmat(v, num_x, 1);

tmp_xu = repmat(xu, num_x, 1);
tmp_xl = repmat(xl, num_x, 1);
tmp_xl(tmp_v) = x;
end