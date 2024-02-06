function[XL, FL, FC, LLcount, ULcount] = LL_dualEvolution(xu, Params, prob, init_pop, vis)
% create initial population

ub = prob.ll_bu;
lb = prob.ll_bl;
LL_nvar = length(ub);

% Expand xu 
funh1 = @(xl)objective_func1(prob, xu, xl);
conh1 = @(xl)constraint_func1(prob, xu, xl);

funh2 = @(xl)objective_func2(prob, xu, xl);
conh2 = @(xl)constraint_func2(prob, xu, xl);

p.gen = Params.LL_gensize;
p.popsize = Params.LL_popsize;

if vis
    pf = prob.PF_LL(129, xu(1, :));
    if isempty(pf)
        pf = [0, 0];
    end

    num = 200;
    llps = prob.PS_LL(num, xu);
    [upf, uc] = prob.evaluate_u(repmat(xu, num, 1), llps);
    id = uc<=0;
    upf = upf(id, :);
else
    pf = [0, 0];
    upf = [0, 0];
end

[XL, FL, FC, archive, ULcount] = ea_dualpopsolver(funh1, funh2, LL_nvar, lb, ub, init_pop, conh1, conh2, p, 'visualize', vis, 'pf', pf, 'upf', upf, 'xu', xu);

% fprintf('LL search takes %d generations \n', max(archive.sols(:, 1))+1)
LLcount = size(archive, 1);

end


%------wrapper----
%----not perfect data structure can cause objective evaluated twice, but it is fine
function  f = objective_func1(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_l(xu, xl);
end

function c = constraint_func1(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_l(xu, xl);
end


function  f = objective_func2(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_u(xu, xl);
end

function c = constraint_func2(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_u(xu, xl);
end

