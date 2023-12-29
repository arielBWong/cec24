function [pop] = reduce_pop(pop, popsize, front_idx, lb, ub)
% Reduction is based on information of solutions that have been sorted
% Note: make sure pop is sorted

select_id = 1 : popsize;  % compatable for SO 

if size(pop.F, 2) > 1      % MO special treatment
    front1_binary = front_idx == 1;
    num_front1 = sum(front1_binary);
    front1_idx = 1: num_front1;   % pop is in sorted status
    
    if num_front1 > popsize
       
        % DSS on F
%         ideal = min(pop.F(front1_binary, :), [], 1);
%         nadir = max(pop.F(front1_binary, :), [], 1);
%         sparseID =  Sparse_Selection_X(front1_idx, pop.F, ideal, nadir, popsize);
%         select_id1 = front1_idx(sparseID);
%         select_id2 = [];
       
        % DSS on X 
        select_id1 = [];
        front1_idx(select_id1) = [];          % sorted pop, 1:num_front1 is the first front
        sparseID = Sparse_Selection(front1_idx, pop.X, lb, ub, popsize);
        select_id2 = front1_idx(sparseID);

        select_id = [select_id1, select_id2];
    end
end

fields = fieldnames(pop);
nf = length(fields);
for i = 1:nf
    if iscell(pop.(fields{i})) % deal with all cell array
        if ~isempty(pop.(fields{i}))
            pop.(fields{i}) = pop.(fields{i})(select_id);
        end
    else
        if ~isempty(pop.(fields{i}))
            pop.(fields{i}) = pop.(fields{i})(select_id, :);
        end
    end
end

return
