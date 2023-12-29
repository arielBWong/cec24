function cx = remove_repeated_solution(existingX, cx)
% cx is child solutions before evaluation
% this method remove new solutions that are already exising in past archive
lia = ismember(cx, existingX, 'rows');
cx(lia, :) = [];
end