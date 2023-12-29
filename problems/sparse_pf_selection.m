function pf = sparse_pf_selection(pf, n)
% this method calls DSS to select close to uniform solutions
% note: to ensure uniformity, n has to be 2^k + 1, k is random integer
candidate_id = 1:size(pf, 1);
LB = min(pf, [], 1);
UB = max(pf, [], 1);
select_id = Sparse_Selection(candidate_id, pf, LB, UB, n);
pf = pf(candidate_id(select_id), :);

end