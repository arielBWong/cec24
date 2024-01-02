function[XL, FL, FC, LLcount] = LL_Evolution(xu, Params, prob, init_pop, vis, termination_condition)
% create initial population

ub = prob.ll_bu;
lb = prob.ll_bl;
LL_nvar = length(ub);

% Expand xu 
funh = @(xl)objective_func(prob, xu, xl);
conh = @(xl)constraint_func(prob, xu, xl);

p.gen = Params.LL_gensize;
p.popsize = Params.LL_popsize;

if vis
    pf = prob.PF_LL(129, xu(1, :));
    if isempty(pf)
        pf = [0, 0];
    end

else
    pf = [0, 0];
end

[XL, FL, FC, archive, ~] = gsolver(funh, LL_nvar, lb, ub, init_pop, conh, p, 'visualize', vis, 'pf', pf, 'termination_criterion', termination_condition);

% fprintf('LL search takes %d generations \n', max(archive.sols(:, 1))+1);
LLcount = size(archive.sols, 1);

end


function plotMO(fighn, pop, pf, gen)
clf(fighn);
plot(pf(:, 1), pf(:, 2), 'go'); hold on;
plot(pop(:, 1), pop(:, 2), 'ro');
title(num2str(gen));
grid on;
pause(0.1);
end

%------wrapper----
%----not perfect data structure can cause objective evaluated twice, but it is fine
function  f = objective_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_l(xu, xl);
end

function c = constraint_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_l(xu, xl);
end