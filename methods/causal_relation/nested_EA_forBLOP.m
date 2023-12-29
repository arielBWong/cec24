function nested_EA_forBLOP(varargin)
% create argument parser
p = inputParser;
addParameter(p, 'prob_str', 'DS5(1, 9)');
addParameter(p, 'strategy', 'net4XLNDinit');
addParameter(p, 'seed', 1);
% parse input parameter 
prob_str = p.Results.prob_str;
strategy = p.Results.strategy;
seed = p.Results.seed;
%-------------------------------


% causal relation determination
prob = eval(prob_str);
[r, xl_causal_fl] = LL_variable_test(prob);

% load parameter 
parameter_setting();
load('parameter.mat');

% run nested EA
XU = unifrnd(repmat(prob.ul_bl, parameter.UL_popsize,1),repmat(prob.ul_bu,parameter.UL_popsize,1));
population = solutions();


for ii = 1: parameter.UL_popsize
    [XL, FL, FLC, ~] = LL_Evolution(XU(ii,:), parameters, prob, [], false, 0); 
    
    % fix XL corresponding UL related variables, if necessary
    if r
        sub_ub = prob.ll_bu(v);
        sub_lb = prob.ll_bl(v);

        for jj = 1:size(XL, 1)
            xl = XL(jj, :);

            p.gen = parameter.LL_gensize;
            p.popsize = parameter.LL_popsize;

            % how to run experiments on UL? use binary v to assign value back
            funh_obj=@(x)partial_search_objective(x, xu, xl, v, prob);
            funh_con=@(x)partial_search_constraint(x, xu, xl, v, prob);
            

            [XL_partial, FU, FC, ~, ~] = gsolver(funh_obj, LL_nvar, sub_lb, sub_ub, [], funh_con, p, 'visualize', true, 'pf', [], 'termination_criterion', 0);

            
        end
    else
        % no action
    end
    
    % continue on UL evolution
    [FU, FC] = prob.evaluate_u(repmat(XU(ii, :), size(XL, 1), 1), XL);
    population.add(XU(ii, :), XL, FU, FL, FC, FLC);
end

% Generate child population 
param_tmp.popsize = Params.UL_popsize;
childXU = generate_child_DE(prob.ul_bl, prob.ul_bu, population.XUs, param_tmp);
childXU = unique(childXU, 'rows', 'stable');
childXU = remove_repeated_solution(unique(Active_pop.XU, 'rows', 'stable'), childXU);






end


function [r, xl_causal_fl] = LL_variable_test(prob)
% causal relation determination
r = false;

xu = unifrnd(prob.ul_bl, prob.ul_bu, 1);
xl = unifrnd(prob.ll_bl, prob.ll_bu, 1);

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


