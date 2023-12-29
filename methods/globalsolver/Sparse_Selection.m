function [selected_ids] = Sparse_Selection(candidate_ids, total_x, LB_x, UB_x, required_size)
% input: 
% candidate_ids: [row] index to pick subset or all solutions from second argument (total_x) for
% selection
% total_x: [2D matrix] targets to be select from
% LB_x: [row] lower bound of total_x (same column size as total_x)
% UB_x: [row] upper bound of total)x (same column size as total)x)
% required_size: how many solutions to be picked up
% output:
% selected_ids: index -> candidate ids for picking up target solutions
% how to use output -> candidate_ids(selected_ids)
%--------
candidate=total_x(candidate_ids,:);

% Normalizing
candidate_norm=(candidate - repmat(LB_x, size(candidate, 1), 1))./(repmat(UB_x-LB_x,size(candidate,1),1));

% This is to detect the first 2 solutions
[D,~] = pdist2(candidate_norm,candidate_norm,'euclidean','Largest',1);
maxD=max(D);

% fighn = figure(1);
% plot(candidate(:, 1), candidate(:, 2), 'k.'); hold on;

% DSS_X
b=find(D==maxD);
idall = [b setdiff(randperm(size(candidate_norm,1)), b)];
selected_ids = idall(1);
leftset = setdiff((1:size(candidate_norm,1)),selected_ids);

while numel(selected_ids) < required_size
    mindist = zeros(numel(leftset),2);
    dist_matrix = pdist2(candidate_norm(leftset, :), candidate_norm(selected_ids, :), 'euclidean');
    dist_matrix = min(dist_matrix, [], 2);
    [~, idmax] = max(dist_matrix);
    selected_ids = [selected_ids leftset(idmax(1))];   
    leftset = setdiff((1:size(candidate_norm,1)),selected_ids);    
    % plot(candidate(selected_ids, 1), candidate(selected_ids, 2), 'ro');
    
end
return
