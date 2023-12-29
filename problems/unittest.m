
prob = DS5m(3,2);
ub = prob.ll_bu;
lb = prob.ll_bl;
LL_nvar = length(ub);

funh = @(xl)objective_func(prob, xu, xl);
conh = @(xl)constraint_func(prob, xu, xl);

xu = [0.5, 1, 1];
p.gen = 300;
p.popsize = 20;
pf = prob.PF_LL(129, xu(1, :));
[bestx, bestf, bestc, archive, external_return] = gsolver(funh, LL_nvar, lb, ub, [], conh, p, 'pf', pf,  'visualize', true, 'termination_criterion', true);



%------wrapper----
%----not perfect data structure can cause objective evaluated twice,
% mitigate by counting evolution solutions
function  f = objective_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[f, ~] = prob.evaluate_l(xu, xl);
end

function c = constraint_func(prob, xu, xl)
xu = repmat(xu, size(xl, 1), 1);
[~, c] = prob.evaluate_l(xu, xl);
end