function [pop, front_idx] = pop_sort(pop)
% auxiliary function
% this function sort evoluation popluation w.r.t.
% number of objectives
% constraints
% nd sort is not compatible with single objective problems
%-----------------------------------------------------
front_idx = [];
numcon = size(pop.C, 2);
numobj = size(pop.F, 2);

if numcon == 0
    % unconstraint problem
    % pop always has member F
    [pop, front_idx] = sort_if_all_feasible(pop, numobj);
    return;
end

%----
if numcon>0
    % process constraints by sum of cv
    cons = pop.C;
    cons(cons <= 0) = 0;
    cons = sum(cons, 2);
    
    infeasible_id = find(cons > 0);
    feasible_id = find(cons <= 0);
    
    infeasible_binary = cons > 0;
    feasible_binary = cons <= 0;
    
    if isempty(infeasible_id) 
        [pop, front_idx] = sort_if_all_feasible(pop, numobj);
    else
        % FF is feasible objectives, CI sumed cv constrait
        FF = pop.F(feasible_binary, :);
        CI = cons(infeasible_binary, :);
        
        if sum(feasible_binary) == 0  %  short cut for no feasible solutions
            [~, ids_FI] = sort(CI);
            infeasible_id = infeasible_id(ids_FI);
            sorted_id = infeasible_id;
            front_idx = ones(size(sorted_id, 1),1) * 2;  % if all infeasible, first one is set front 1
            front_idx(1) = 1;
        else
            % process feasible
            if numobj > 1                                               % mo problem
                [fronts, ids_FF, ~] = nd_sort(FF, (1:size(FF, 1))');
                for j = 1:size(fronts,2)
                    front_idx = [front_idx; ones(length(fronts(j).f), 1) * j];
                end
            else                                                        % so problem
                [~, ids_FF] = sort(FF);                                 % acending sort/minimization
                front_idx = ones(length(ids_FF), 1);
            end
            feasible_id = feasible_id(ids_FF);
            
            % process infeasible
            [~, ids_FI] = sort(CI);
            infeasible_id = infeasible_id(ids_FI);
            num_front = max(front_idx);
            front_idx = [front_idx; ones(length(ids_FI), 1) * (num_front+1)];
            
            %-----combine id----
            sorted_id = [feasible_id; infeasible_id];    
        end
            
        %--------------
        fields = fieldnames(pop);
        for i = 1:length(fields)
            content = pop.(fields{i});
            if ~isempty(content)
                pop.(fields{i}) = content(sorted_id,:); % adjust order
            end
        end  
    end
    
end 

end

function [pop, front_idx] = sort_if_all_feasible(pop, numobj)
% common method for sorting
front_idx = []; % indicating which front, the current solution is in
if numobj > 1                                               % mo problem
    [fronts, ids, ~] = nd_sort(pop.F, (1:size(pop.F, 1))');
    
    for j = 1:size(fronts,2)
        front_idx = [front_idx; ones(length(fronts(j).f), 1)*j];
    end
    
else                                                        % so problem
    [~, ids] = sort(pop.F);                                 % acending sort/minimization
    front_idx = ones(size(pop.F, 1), 1) * 2;                % So only best solution is front 1, rest are front2
    front_idx(1) = 1;
end

fields = fieldnames(pop);
for i = 1:length(fields)
    content = pop.(fields{i});
    if ~isempty(content)
        pop.(fields{i}) = content(ids,:); % adjust order
    end
end
end

